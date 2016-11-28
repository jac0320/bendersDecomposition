/*
 * master.c
 *
 *  Created on: Nov 28, 2016
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University.
 *  	Report: harsha(at)smu(dot)edu
 */

#include "benders.h"
configType config;

int solveMaster(probType *prob, cellType *cell) {
	int 	status, n;

	if ( config.MASTER_TYPE == PROB_QP ) {
		/* In regularized QP method, switch to barrier optimizer to solve quadratic master problem */
		changeSolverType(ALG_BARRIER);

#ifdef MODIFY_CHECK
		writeProblem(cell->master->lp, "master_k.lp");
#endif

		if ( solveProblem(cell->master->lp, cell->master->name, PROB_QP, &status) ) {
			if ( status == STAT_INFEASIBLE ) {
				/*  master is infeasible, terminate the algorithm */
				printf("Master is infeasible, terminating the algorithm.\n");
				writeProblem(cell->subprob->lp, "infeasMaster.lp");
				return 1;
			}
			else {
				errMsg("algorithm", "solveMaster", "failed to solve problem in solver", 0);
				return 1;
			}
		}

		if ( getPrimal(cell->master->lp, cell->candidU, prob->num->prevCols) ) {
			errMsg("algorithm", "solveMaster", "failed to obtain the optimal primal solution", 0);
			return 1;
		}

		cell->candidU[0] = 0.0;
		for ( n = 1; n <= prob->num->prevCols; n++ )
			cell->candidU[n] += cell->incumbU[n];
		cell->candidU[0] = oneNorm(cell->candidU, prob->num->prevCols);

		/* Switch back to primal/dual simplex solver */
		changeSolverType(ALG_AUTOMATIC);

		/* Note improvement at candidate solution with respect to current approximation */
		cell->improve = maxCutHeight(cell->cuts, cell->candidU, prob->coord->colsC, prob->num->cntCcols) - cell->incumbEst;
	}

	return 0;
}//END solveMaster()

void checkImprovement(probType **prob, cellType *cell) {
	double candidEst;

	cell->incumbEst = vXv(cell->master->objx, cell->incumbU, NULL, prob[0]->num->cols) +
			maxCutHeight(cell->cuts, cell->incumbU, prob[1]->coord->colsC, prob[1]->num->cntCcols);
	candidEst = vXv(cell->master->objx, cell->candidU, NULL, prob[0]->num->cols) +
			maxCutHeight(cell->cuts, cell->candidU, prob[1]->coord->colsC, prob[1]->num->cntCcols);

#ifdef ALGO_TRACE
	printf("Incumbent estimate = %lf\tCandidate estimate = %lf\n", cell->incumbEst, candidEst);
#endif

	if ( (candidEst - cell->incumbEst) < config.R1*cell->improve) {
		/* incumbent change is recommended */
		if ( replaceIncumbent(prob, cell, candidEst) )
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent", 0);
		printf("+"); fflush(stdout);
	}

}//END checkImprovement()

int replaceIncumbent(probType **prob, cellType *cell, double candidEst) {
	int n;

	/* replace the incumbent solution with the current candidate solution */
	copyVector(cell->candidU, cell->incumbU, prob[0]->num->cols, 1);

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( changeQPrhs(cell->master->lp, prob[1]->coord->colsC, prob[1]->num->cntCcols, prob[0]->num->rows, prob[0]->Dbar, prob[0]->bBar, cell->cuts, cell->incumbU) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}
	if ( changeQPbds(cell->master->lp, prob[0]->num->cols, prob[0]->sp->bdl, prob[0]->sp->bdu, cell->incumbU) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the bounds after incumbent update", 0);
		return 1;
	}



	cell->incumbEst = candidEst;
	cell->improve	= 0.0;

	return 0;
}//END replaceIncumbent()

int constructQP(LPptr lp, int numCols, double regSigma) {
	vector	qpsepvec;
	int		n, status;

	/* allocate memory for regularizing terms */
	if ( !(qpsepvec = (vector) arr_alloc(numCols+1, double)) )
		errMsg("allocation", "constructQP", "regularization term vector", 0);

	/* (i) Add regularization term */
	/* matrix with diagonal elements set to $\sigma/2$, the last term corresponds to $\eta$ */
	for ( n = 0; n < numCols; n++ )
		qpsepvec[n] = regSigma;
	qpsepvec[n] = 0;

	/* Now copy the Q matrix for in the solver. */
	status = copyQPseparable(lp, qpsepvec);
	if ( status ) {
		errMsg("solver", "constructQP", "copy Q matrix into the solver", 0);
		return 1;
	}

	mem_free(qpsepvec);
	return 0;
}//END constructQP()

int changeQPrhs(LPptr lp, intvec betaCols, int betaLen, int numRows, sparseMatrix *Dbar, sparseVector *bBar, cutsType *cuts, vector X) {
	vector 	qpRHS;
	intvec indices;
	int		n;

	if ( !(qpRHS = (vector) arr_alloc(numRows+1, double)))
			errMsg("allocation", "computeRHS", "indices", 0);
	if ( !(indices = (intvec) arr_alloc(numRows, int)))
		errMsg("allocation", "computeRHS", "indices", 0);

	/* change the right-hand side to \bar{b}_t - D_t \hat{u}_t */
	for ( n = 0; n < numRows; n++)
		indices[n]= n;

	/* for root stage copy the right-hand side fixed part */
	for ( n = 1; n <= bBar->cnt; n++ )
		qpRHS[bBar->col[n]] += bBar->val[n];

	qpRHS = MSparsexvSub(Dbar, X, qpRHS);

	/* changing cut right-hand side of stage subproblem */
	for ( n = 0; n < cuts->cnt; n++ )
		qpRHS[cuts->vals[n]->rowNum+1] = cuts->vals[n]->alpha - vXv(cuts->vals[n]->beta, X, betaCols, betaLen);

	mem_free(indices); mem_free(qpRHS);

	return 0;
}//END changeQPrhs()

int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector X) {
	intvec 	Cindices;
	vector	lbounds, ubounds;
	string	llu, ulu;
	int		n;

	if (!(Cindices = (intvec) arr_alloc(numCols, int)))
		errMsg ("allocation", "change_bounds", "Cindices", 0);
	if (!(lbounds = (vector) arr_alloc(numCols, double)))
		errMsg ("allocation", "change_bounds", "lbounds", 0);
	if (!(llu= (string) arr_alloc(numCols, char)))
		errMsg ("allocation", "change_bounds", "llu", 0);

	if (!(ubounds = (vector) arr_alloc(numCols, double)))
		errMsg ("allocation", "change_bounds", "ubounds", 0);
	if (!(ulu= (string) arr_alloc(numCols, char)))
		errMsg ("allocation", "change_bounds", "ulu", 0);

	for (n = 0; n < numCols; n++) {
		Cindices[n] = n;
		ubounds[n] = bdu[n]- X[n+1];
		ulu[n] = 'U';
		lbounds[n] = bdl[n] - X[n+1];
		llu[n] = 'L';
	}

	if ( changeBDS (lp, numCols, Cindices, ulu, ubounds) ) {
		errMsg("solver", "updtQPbds", "failed to change the upper bounds with incumbent information", 0);
		return 1;
	}

	if ( changeBDS (lp, numCols, Cindices, llu, lbounds) ) {
		errMsg("solver", "updtQPbds", "failed to change the lower bounds with incumbent information", 0);
		return 1;
	}

	mem_free(Cindices);
	mem_free(lbounds); mem_free(llu);
	mem_free(ubounds); mem_free(ulu);

	return 0;
}//END changeQPbds()

double maxCutHeight(cutsType *cuts, vector x, intvec betaIdx, int betaLen) {
	double 	val, ht = -DBL_MAX;
	int		n;

	for ( n = 0; n < cuts->cnt; n++ ) {
		val = cuts->vals[n]->alpha - vXv(cuts->vals[n]->beta, x, betaIdx, betaLen);
		if ( ht < val )
			ht = val;
	}

	return val;
}//END maxCutHeight()

oneProblem *newMaster(probType *prob, cutsType *cuts, omegaType *omega) {
	oneProblem *master = NULL;
	char		*q, cutName[NAMESIZE] = { "Old    " }, etaName[NAMESIZE];
	int			j, cnt, idx, rowOffset, colOffset, cutsCnt = 0;

	if ( cuts != NULL )
		cutsCnt = cuts->cnt;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("allocation", "newMaster", "copy", 0);
	master->bdl  = NULL; master->bdu 	 = NULL; master->cname  = NULL; master->cstore = NULL; master->ctype  = NULL;
	master->lp   = NULL; master->matbeg  = NULL; master->matcnt = NULL; master->matind = NULL; master->matval = NULL;
	master->name = NULL; master->objname = NULL; master->objx   = NULL; master->rhsx   = NULL; master->rname  = NULL;
	master->rstore=NULL; master->senx    = NULL;

	/* assign values to scalar terms */
	master->objsen = prob->sp->objsen;
	master->marsz = prob->sp->mar + cutsCnt;	master->mar = prob->sp->mar + cutsCnt;	master->mac = prob->sp->mac; /* Note: mac does not include eta columms */
	if ( config.MULTICUT ) {
		/* Additional columns to hold eta's (surrogates for individual recourse functions), one for each possible outcome */
		master->macsz = prob->sp->mac + omega->cnt; master->cstorsz = prob->sp->cstorsz + omega->cnt*NAMESIZE;
	}
	else {
		/* one additional column to hold eta (surrogate for expected recourse function */
		master->macsz = prob->sp->mac+1; master->cstorsz = prob->sp->cstorsz + NAMESIZE;
	}
	master->matsz  = prob->sp->matsz + (cutsCnt)*(master->mac + 1);
	master->numnz  = prob->sp->numnz;
	master->rstorsz = prob->sp->rstorsz + NAMESIZE*(cutsCnt);
	master->numInt = prob->sp->numInt; master->numBin = prob->sp->numBin;
	if ( master->numInt + master->numBin > 0 ) {
		if ( config.PROXIMAL == 0) {
			if ( config.MASTER_TYPE == PROB_LP || config.MASTER_TYPE == PROB_QP)
				printf("Warning :: Ignoring integrality and solving a relaxed problem.\n");
			master->type = config.MASTER_TYPE;
		}
		else {
			if ( config.MASTER_TYPE == PROB_MILP || config.MASTER_TYPE == PROB_MIQP ) {
				errMsg("configuration", "newMaster", "No regularization supported for MILPs and MIQPs", 0);
				return NULL;
			}
			printf("Warning :: Ignoring integrality and solving a regularized, relaxed problem.\n");
			master->type = PROB_QP;
		}
	}
	else
		(config.PROXIMAL == 1) ? (master->type = PROB_QP): (master->type = config.MASTER_TYPE);

	/* Make all allocations of known sizes, as calculated above */
	if (!(master->name = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newMaster", "master->name",0);
	if (!(master->objname = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newMaster", "master->objname",0);
	if (!(master->objx = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->objx",0);
	if (!(master->bdl = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->bdl",0);
	if (!(master->bdu = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->bdu",0);
	if (!(master->rhsx = arr_alloc(master->marsz, double)))
		errMsg("allocation", "newMaster", "master->rhsx",0);
	if (!(master->senx = arr_alloc(master->marsz, char)))
		errMsg("allocation", "newMaster", "master->senx",0);
	if (!(master->matbeg = arr_alloc(master->macsz, int)))
		errMsg("allocation", "newMaster", "master->matbeg",0);
	if (!(master->matcnt = arr_alloc(master->macsz, int)))
		errMsg("allocation", "newMaster", "master->matcnt",0);
	if (!(master->rname = arr_alloc(master->marsz, string)))
		errMsg("allocation", "newMaster", "master->rname",0);
	if (!(master->rstore = arr_alloc(master->rstorsz, char)))
		errMsg("allocation", "newMaster", "master->rstore",0);
	if (!(master->cname = arr_alloc(master->macsz, string)))
		errMsg("allocation", "newMaster", "master->cname",0);
	if (!(master->cstore = arr_alloc(master->cstorsz, char)))
		errMsg("allocation", "newMaster", "master->cstore",0);
	if (!(master->matval = arr_alloc(master->matsz, double)))
		errMsg("allocation", "newMaster", "master->matval",0);
	if (!(master->matind = arr_alloc(master->matsz, int)))
		errMsg("allocation", "newMaster", "master->matind",0);
	if (!(master->ctype = arr_alloc(master->macsz, char)))
		errMsg("allocation", "newMaster", "master->ctype",0);

	/* First copy information directly from the original subproblem. */
	/* Copy the subproblem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	j = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		master->cstore[j++] = *q;

	j = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		master->rstore[j++] = *q;

	strcpy(master->name, prob->sp->name);
	strcpy(master->objname, prob->sp->objname);

	/* Calculate difference in pointers for oneProblems in probType and cellType for row and column names */
	colOffset = master->cstore - prob->sp->cname[0];
	rowOffset = master->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original subproblem */
	master->numnz = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		master->objx[j] = prob->sp->objx[j];
		master->ctype[j] = prob->sp->ctype[j];
		master->bdu[j] = prob->sp->bdu[j];
		master->bdl[j] = prob->sp->bdl[j];
		master->cname[j] = prob->sp->cname[j] + colOffset;
		master->matbeg[j] = master->numnz;
		master->matcnt[j] = prob->sp->matcnt[j];
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			master->matval[master->numnz] = prob->sp->matval[idx];
			master->matind[master->numnz] = prob->sp->matind[idx];
			master->numnz++;
		}
	}

	/* Copy all row information from the original subproblem */
	for (j = 0; j < prob->sp->mar; j++) {
		master->rhsx[j] = prob->sp->rhsx[j];
		master->senx[j] = prob->sp->senx[j];
		master->rname[j] = prob->sp->rname[j] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	if (config.MULTICUT == 0) {
		strcpy(master->cstore + prob->sp->cstorsz, "eta");
		master->cname[master->mac] = master->cstore + prob->sp->cstorsz;
		master->objx[master->mac] = 1.0;
		master->ctype[master->mac] = 'C';
		master->bdu[master->mac] = INFBOUND;
		master->bdl[master->mac] = prob->lb;
		master->matbeg[master->mac] = master->numnz;
		master->matcnt[master->mac] = cutsCnt;
	}
	else {
		colOffset = prob->sp->cstorsz;
		for (cnt = 0; cnt < omega->cnt; cnt ++) {
			sprintf(etaName, "eta%d", cnt + 1);
			strcat(master->cstore + colOffset, etaName);
			master->cname[master->mac + cnt] = master->cstore + colOffset;
			master->objx[master->mac + cnt] = omega->probs[cnt];
			master->ctype[master->mac + cnt] = 'C';
			master->bdu[master->mac + cnt] = INFBOUND;
			master->bdl[master->mac + cnt] = prob->lb;
			master->matbeg[master->mac + cnt] = master->numnz;
			master->matcnt[master->mac + cnt] = 0;
			colOffset += strlen(etaName)+1;
		}
	}

	if (cutsCnt > 0) {
		/* Now copy information from the cuts into the new master problem. */
		for (idx = 0; idx < cuts->cnt; idx++) {
			master->matval[idx + cnt] = 1.0;
			master->matind[idx + cnt] = master->mar + idx;
		}
		/* Copy the constraint coefficients */
		for (j = 0; j < master->mac; j++) {
			cnt = master->matbeg[j] + master->matcnt[j];
			master->matcnt[j] += cuts->cnt;

			for (idx = 0; idx < cuts->cnt; idx++) {
				master->matval[cnt + idx] = cuts->vals[idx]->beta[j + 1];
				master->matind[cnt + idx] = master->mar + idx;
			}
		}

		/* Give names to the cut constraints, and add rhs values */
		idx = master->rstorsz;
		for (cnt = 0; cnt < cuts->cnt; cnt++) {
			cutName[3] = '0' + cnt / 1000 % 10;
			cutName[4] = '0' + cnt / 100 % 10;
			cutName[5] = '0' + cnt / 10 % 10;
			cutName[6] = '0' + cnt / 1 % 10;
			master->rname[master->mar + cnt] = master->rstore + idx;
			strcpy(master->rname[master->mar + cnt], cutName);
			idx += strlen(cutName) + 1;
			master->rhsx[master->mar + cnt] = cuts->vals[cnt]->alpha;
			master->senx[master->mar + cnt] = 'G';
		}
	}

	/* Load the copy into CPLEX now */
	master->lp = setupProblem(master->name, master->type, master->macsz, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,
			master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "newMaster", "failed to setup master problem in the solver",0);
		return NULL;
	}

#ifdef INPUT_CHECK
	writeProblem(master->lp, "master.lp");
#endif

	return master;
}//END newMaster()
