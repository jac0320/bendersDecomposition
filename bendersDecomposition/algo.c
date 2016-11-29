/*
 * algo.c
 *
 *  Created on: Nov 22, 2016
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *      e-mail: harsha(at)smu(dot)edu 
 */


#include "benders.h"

configType config;

int algo(oneProblem *orig, stocType *stoc, timeType *tim) {
	probType **prob;
	cellType *cell;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) ) {
		errMsg("algorithm", "algo", "failed to setup the algorithm",0);
		goto TERMINATE;
	}

	/* main loop of the algorithm */
	while (cell->k < config.MAX_ITER) {
#ifdef ALGO_TRACE
		printf("\nIterarion-%3d ::", cell->k+1);
#else
		if (cell->k % 100 == 0)
			printf("\nIterarion-%3d ::", cell->k+1);
#endif
		cell->k++;

		/* Step 1: Solve the subproblems and create the optimality cut*/
		if ( solveSubprobs(prob[1], cell) ) {
			errMsg("algorithm", "algo", "failed to solve the subproblem",0);
			goto TERMINATE;
		}

		/* Step 2: Cut generation */
		if ( config.MULTICUT == 0 ) {
			if ( formSingleCut(prob, cell) ) {
				errMsg("algorithm", "algo", "failed to generate a new cut",0);
				goto TERMINATE;
			}
		}
		else {
			if ( formMultiCut(prob, cell) ) {
				errMsg("algorithm", "algo", "failed to generate multi cuts",0);
				goto TERMINATE;
			}
		}

		/* Step 3: Check improvement to determine if an incumbent update is necessary */
		if (config.PROXIMAL == 1 )
			checkImprovement(prob, cell);

		/* Step 4: Solve the master problem */
		if ( solveMaster(prob[1], cell) ) {
			errMsg("algorithm", "algo", "failed to solve the master problem",0);
			goto TERMINATE;
		}

	}//END main loop

	freeCellType(cell);
	freeProbType(prob, tim->numRows);
	return 0;

	TERMINATE:
	freeCellType(cell);
	freeProbType(prob, tim->numRows);
	return 1;
}//END algo()

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell) {
	vector	u0 = NULL, lb = NULL;

	/* setup mean value problem which will act as reference for all future computations */
	u0 = meanProblem(orig, stoc);
	if ( u0 == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		goto TERMINATE;
	}

	/* obtain a lower bound on optimal objective function value */
	lb = calcLowerBound(orig, tim);
	if ( lb == NULL ) {
		errMsg("setup", "setupAlgo", "failed to obtain the lower bound", 0);
		goto TERMINATE;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.MAX_ITER);
	if ((*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to decompose the problem into stage problems", 0);
		goto TERMINATE;
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell((*prob), stoc, u0);
	if ((*cell) == NULL) {
		errMsg("setup", "setupAlgo", "failed to setup the cell used by the algorithm", 0);
		goto TERMINATE;
	}

	mem_free(lb);
	mem_free(u0);
	return 0;

	TERMINATE:
	if (lb) mem_free(lb);
	if (u0) mem_free(u0);
	return 1;
}//END setupAlgo()

cellType *newCell(probType **prob, stocType *stoc, vector u0) {
	cellType *cell;

	if (!(cell = (cellType *) mem_malloc (sizeof(cellType))))
		errMsg("allocation", "newCell", "cell",0);
	cell->k = 0;

	/* stochastic elements */
	cell->omega   = newOmega(stoc);
	cell->sigma   = newSigma(config.MAX_ITER, cell->omega->cnt);
	cell->delta   = newDelta(config.MAX_ITER, cell->omega->cnt);

	/* oneProblems for master and subproblem */
	cell->master  = newMaster(prob[0], NULL, cell->omega);
	cell->subprob = newSubproblem(prob[1]);

	cell->candidU = duplicVector(u0, prob[1]->num->cols);

	if (config.PROXIMAL == 1) {
		cell->incumbU = duplicVector(u0, prob[1]->num->cols);
		cell->incumbEst = vXv(cell->master->objx, cell->incumbU, NULL, prob[0]->num->cols);
		if ( config.MULTICUT == 1) {
			cell->maxCuts = config.CUT_MULT*prob[0]->num->cols + cell->omega->cnt + 1;
		}
		else {
			cell->maxCuts = config.CUT_MULT*prob[0]->num->cols + 1;
		}
		if ( !(cell->masterPi = (vector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
			errMsg("allocation", "newCell", "cell->masterPi", 0);
		cell->improve = 0.0;
		cell->quadScalar = config.QUAD_SCALAR;
	}
	else {
		cell->incumbU = NULL;
		cell->incumbEst = 0.0;
		cell->maxCuts = config.MAX_ITER;
		cell->masterPi = NULL;
	}

	if ( !(cell->cuts = (cutsType *) mem_malloc(sizeof(cutsType))) )
		errMsg("allocation", "newCell", "cell->cuts", 0);
	if ( !(cell->cuts->vals = (oneCut **) arr_alloc(cell->maxCuts, oneCut *)) )
		errMsg("allocation", "newCell", "cell->cuts->vals", 0);
	cell->cuts->cnt = 0;

	if ( config.PROXIMAL == 1 ) {
		/* construct the proximal term */
		if ( constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to setup the proximal term", 0);
			return NULL;
		}

		/* update the right-hand side and the bounds with new incumbent solution */
		if ( changeQPrhs(cell->master->lp, prob[1]->coord->colsC, prob[1]->num->cntCcols, prob[0]->num->rows, prob[0]->Dbar, prob[0]->bBar, cell->cuts, cell->incumbU) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}
		if ( changeQPbds(cell->master->lp, prob[0]->num->cols, prob[0]->sp->bdl, prob[0]->sp->bdu, cell->incumbU) ) {
			errMsg("setup", "newCell", "failed to change the bounds after incumbent update", 0);
			return NULL;
		}

	}

	return cell;
}//END newCell()

oneProblem *newSubproblem(probType *prob) {
	oneProblem  *subprob = NULL;
	char		*q;
	int			j, cnt, idx, rowOffset, colOffset;

	if (!(subprob = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("allocation", "newSubproblem", "copy", 0);
	subprob->bdl  = NULL; subprob->bdu 	  = NULL; subprob->cname  = NULL; subprob->cstore = NULL; subprob->ctype = NULL;
	subprob->lp   = NULL; subprob->matbeg = NULL; subprob->matcnt = NULL; subprob->matind = NULL; subprob->matval = NULL;
	subprob->name = NULL; subprob->objname= NULL; subprob->objx   = NULL; subprob->rhsx	  = NULL; subprob->rname = NULL;
	subprob->rstore=NULL; subprob->senx = NULL;

	/* assign values to scalar terms */
	subprob->objsen = prob->sp->objsen;
	subprob->matsz = prob->sp->matsz;		subprob->numnz = prob->sp->numnz;
	subprob->marsz = prob->sp->mar;			subprob->mar = prob->sp->mar;
	subprob->macsz = prob->sp->mac; 		subprob->mac = prob->sp->mac;
	subprob->cstorsz = prob->sp->cstorsz;	subprob->rstorsz = prob->sp->rstorsz;
	subprob->numInt = prob->sp->numInt;		subprob->numBin = prob->sp->numBin;
	if ( subprob->numInt + subprob->numBin > 0 ) {
		subprob->type = PROB_MILP;
		errMsg("read", "newSubproblem", "subproblem has discrete variables", 0);
		return NULL;
	}
	else
		subprob->type = PROB_LP;

	/* Make all allocations of known sizes, as calculated above */
	if (!(subprob->name = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newSubproblem", "subprob->name",0);
	if (!(subprob->objname = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newSubproblem", "subprob->objname",0);
	if (!(subprob->objx = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->objx",0);
	if (!(subprob->bdl = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->bdl",0);
	if (!(subprob->bdu = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->bdu",0);
	if (!(subprob->rhsx = arr_alloc(subprob->marsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->rhsx",0);
	if (!(subprob->senx = arr_alloc(subprob->marsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->senx",0);
	if (!(subprob->matbeg = arr_alloc(subprob->macsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matbeg",0);
	if (!(subprob->matcnt = arr_alloc(subprob->macsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matcnt",0);
	if (!(subprob->cname = arr_alloc(subprob->macsz, string)))
		errMsg("allocation", "newSubproblem", "subprob->cname",0);
	if (!(subprob->cstore = arr_alloc(subprob->cstorsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->cstore",0);
	if (!(subprob->rname = arr_alloc(subprob->marsz, string)))
		errMsg("allocation", "newSubproblem", "subprob->rname",0);
	if (!(subprob->rstore = arr_alloc(subprob->rstorsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->rstore",0);
	if (!(subprob->matval = arr_alloc(subprob->matsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->matval",0);
	if (!(subprob->matind = arr_alloc(subprob->matsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matind",0);
	if (!(subprob->ctype = arr_alloc(subprob->macsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->ctype",0);

	/* First copy information directly from the original subproblem. */
	/* Copy the subproblem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	j = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		subprob->cstore[j++] = *q;

	j = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		subprob->rstore[j++] = *q;

	strcpy(subprob->name, prob->sp->name);
	strcpy(subprob->objname, prob->sp->objname);

	/* Calculate difference in pointers for oneProblems in probType and cellType for row and column names */
	colOffset = subprob->cstore - prob->sp->cname[0];
	rowOffset = subprob->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original subproblem */
	cnt = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		subprob->objx[j] = prob->sp->objx[j];
		subprob->ctype[j] = prob->sp->ctype[j];
		subprob->bdu[j] = prob->sp->bdu[j];
		subprob->bdl[j] = prob->sp->bdl[j];
		subprob->cname[j] = prob->sp->cname[j] + colOffset;
		subprob->matbeg[j] = cnt;
		subprob->matcnt[j] = prob->sp->matcnt[j];
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			subprob->matval[cnt] = prob->sp->matval[idx];
			subprob->matind[cnt] = prob->sp->matind[idx];
			cnt++;
		}
	}

	/* Copy all row information from the original subproblem */
	for (j = 0; j < prob->sp->mar; j++) {
		subprob->rhsx[j] = prob->sp->rhsx[j];
		subprob->senx[j] = prob->sp->senx[j];
		subprob->rname[j] = prob->sp->rname[j] + rowOffset;
	}

	/* load the information on to the solver */
	subprob->lp = setupProblem(subprob->name, subprob->type, subprob->mac, subprob->mar, subprob->objsen, subprob->objx, subprob->rhsx, subprob->senx,
			subprob->matbeg, subprob->matcnt, subprob->matind, subprob->matval, subprob->bdl, subprob->bdu, NULL, subprob->cname, subprob->rname,
			subprob->ctype);
	if ( subprob->lp == NULL ) {
		errMsg("solver", "newSubprob", "subprob",0);
		return NULL;
	}

#ifdef INPUT_CHECK
	writeProblem(subprob->lp, "subproblem.lp");
#endif

	return subprob;
}//END newSubproblem()

omegaType *newOmega(stocType *stoc) {
	omegaType 	*omega;
	double		val, cumm;
	int 		cnt, i, j, maxObs;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation","newOmega", "omega", 0);

	/* check to see if the problem can be solved as with true distribution */
	// TODO: More detailed check should be included - continuous distribution, indep with small sized scenario set, etc.
	if ( config.SAA == 0 && strstr(stoc->type, "INDEP") != NULL ) {
		cnt = 0; maxObs = 1;
		while ( cnt < stoc->numOmega && maxObs <= 10000) {
			maxObs *= stoc->numVals[cnt];
			cnt++;
		}
	}
	if ( config.SAA == 1 || maxObs > 10000 ) {
		/* SAA */
		printf("Using SAA configuration or maximum number of possible observations exceed 1 million.\nEnter the number of observations to be used in SAA : ");
		scanf("%d", &config.SAA_OBS);

		if ( !(omega->probs = (vector) arr_alloc(config.SAA_OBS, double)))
			errMsg("allocation", "newOmega", "omega->probs", 0);
		if ( !(omega->vals = (vector *) arr_alloc(config.SAA_OBS, vector)) )
			errMsg("allocation", "newOmega", "omega->vals", 0);
		if ( !(omega->istar = (intvec) arr_alloc(config.SAA_OBS, int)))
			errMsg("allocation", "newOmega", "omega->istar", 0);
		omega->numRV = stoc->numOmega;
		omega->cnt = config.SAA_OBS;
		for ( cnt = 0; cnt < omega->cnt; cnt++) {
			omega->probs[cnt]= (1/(double) omega->cnt);
			if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
				errMsg("allocation", "newOmega", "omega->vals[cnt]", 0);
			for (i = 0; i < omega->numRV; i++) {
				val = scalit(0, 1, &config.RUN_SEED); cumm = 0.0;
				for (j = 0; val > cumm; j++)
					cumm += stoc->probs[i][j];
				omega->vals[cnt][i+1] = stoc->vals[i][j] - stoc->mean[i];
			}
		}
	}
	else if ( strstr(stoc->type, "INDEP") != NULL ) {
		if ( !(omega->probs = (vector) arr_alloc(maxObs, double)) )
			errMsg("allocation", "newOmega", "omega->probs", 0);
		if ( !(omega->vals = (vector *) arr_alloc(maxObs, vector)) )
			errMsg("allocation", "newOmega", "omega->vals", 0);
		if ( !(omega->istar = (intvec) arr_alloc(maxObs, int)))
			errMsg("allocation", "newOmega", "omega->istar", 0);
		omega->numRV = stoc->numOmega;
		omega->cnt = maxObs;
		/* all possible observations are generated by permutation with repetition. */
		for ( cnt = 0; cnt < omega->cnt; cnt++) {
			if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
				errMsg("allocation", "newOmega", "omega->vals[n]", 0);
			j = 1; omega->probs[cnt] = 1.0;
			for ( i = omega->numRV-1; i >= 0; i-- ) {
				omega->vals[cnt][i+1] = stoc->vals[i][cnt/(j) % stoc->numVals[i]];
				omega->probs[cnt] *= stoc->probs[i][cnt/(j) % stoc->numVals[i]];
				j *= stoc->numVals[i];
			}
		}
	}
	else if ( strstr(stoc->type, "BLOCKS") != NULL ) {
		/* BLOCKS */
		if ( !(omega->probs = (vector) arr_alloc(stoc->numVals[0], double)))
			errMsg("allocation", "newOmega", "omega->probs", 0);
		if ( !(omega->vals = (vector *) arr_alloc(stoc->numVals[0], vector)) )
			errMsg("allocation", "newOmega", "omega->vals", 0);
		if ( !(omega->istar = (intvec) arr_alloc(stoc->numVals[0], int)))
			errMsg("allocation", "newOmega", "omega->istar", 0);
		omega->numRV = stoc->numOmega;
		omega->cnt=stoc->numVals[0];
		for ( cnt = 0; cnt < omega->cnt; cnt++) {
			omega->probs[cnt]= stoc->probs[0][cnt];
			if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
				errMsg("allocation", "newOmega", "omega->vals[cnt]", 0);
			for (i = 0; i < omega->numRV; i++)
				omega->vals[cnt][i+1] = stoc->vals[i][cnt] - stoc->mean[i];
			omega->vals[cnt][0] = oneNorm(omega->vals[cnt]+1, omega->numRV);
		}
	}
	else {
		errMsg("setup", "newOmega", "setup for other stoch types is under development", 0);
		return NULL;
	}

	return omega;
}//END newOmega()

sigmaType *newSigma(int numIter, int maxObs) {
	sigmaType 	*sigma;

	if (!(sigma = (sigmaType *) mem_malloc(sizeof(sigmaType))) )
		errMsg("allocation", "newSigma", "sigma structure", 0);
	if ( !(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newSigma", "sigma->lambdaIdx", 0);
	if ( !(sigma->vals = (pixbCType *) arr_alloc(numIter, pixbCType)) )
		errMsg("allocation", "newSigma", "sigma->vals", 0);
	sigma->cnt = 0;

	return sigma;
}//END newSigma()

deltaType *newDelta(int numIter, int maxObs) {
	deltaType *delta;

	if (!(delta = (deltaType *) mem_malloc(sizeof(deltaType))) )
		errMsg("allocation", "newDelta", "delta", 0);
	if ( !(delta->lambda = (vector *) arr_alloc(maxObs*numIter, vector)) )
		errMsg("allocation", "newLambda", "delta->vals", 0);
	if ( !(delta->vals = (pixbCType **) arr_alloc(maxObs*numIter, pixbCType *)) )
		errMsg("allocation", "newDelta", "delta->vals", 0);
	delta->rows = 0;
	delta->cols = maxObs;

	return delta;
}//END newDelta()

void freeCellType(cellType *cell) {

	if (cell) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->subprob) freeOneProblem(cell->subprob);
		if (cell->delta) freeDeltaType(cell->delta, cell->omega->cnt);
		if (cell->sigma) freeSigmaType(cell->sigma);
		if (cell->omega) freeOmegaType(cell->omega);
		if (cell->candidU) mem_free(cell->candidU);
		if (cell->incumbU) mem_free(cell->incumbU);
		if (cell->masterPi) mem_free(cell->masterPi);
		if (cell->cuts) freeCutsType(cell->cuts);
		mem_free(cell);
	}

}//END freeCellType()

void freeOmegaType(omegaType *omega) {
	int n;

	if (omega->probs) mem_free(omega->probs);
	if (omega->vals) {
		for ( n = 0; n < omega->cnt; n++ )
			if (omega->vals[n]) mem_free(omega->vals[n]);
		mem_free(omega->vals);
	}
	if (omega->istar) mem_free(omega->istar);
	mem_free(omega);

}//END freeOmegaType()

void freeSigmaType(sigmaType *sigma) {
	int n;

	if (sigma) {
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		if (sigma->vals) mem_free(sigma->vals);
		mem_free(sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int omegaCnt) {
	int m, n;

	if (delta) {
		if (delta->vals) {
			for ( m = 0; m < delta->rows; m++ ) {
				if (delta->vals[m])
//					for ( n = 0; n < omegaCnt; n++ )
//						mem_free(delta->vals[m][n].piC);
					mem_free(delta->vals[m]);
			}
			mem_free(delta->vals);
		}
		if (delta->lambda){
			for ( n = 0; n < delta->rows; n++ ) {
				if (delta->lambda[n])
					mem_free(delta->lambda[n]);
			}
			mem_free(delta->lambda);
		}
		mem_free(delta);
	}

}//END freeDeltaType()
