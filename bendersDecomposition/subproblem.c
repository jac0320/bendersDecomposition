/*
 * subproblem.c
 *
 *  Created on: Nov 24, 2016
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *      e-mail: harsha(at)smu(dot)edu 
 */

#include "benders.h"

configType config;

int solveSubprobs(probType *prob, cellType *cell) {
	vector	rhsx;
	int 	obs, cnt, status;
#ifdef ALGO_TRACE
	double  val, cummVal = 0.0;
	int 	idxSigma, idxLambda;
#endif

	if ( !(rhsx = (vector) arr_alloc(prob->num->rows+1, double)) )
		errMsg("allocation", "solveSubprobs", "rhsx", 0);

	/* 1. copy the original right-hand side and update it with endogenous information (master solution) */
	for (cnt = 1; cnt <= prob->bBar->cnt; cnt++)
		rhsx[prob->bBar->col[cnt]] = prob->bBar->val[cnt];
	rhsx = MSparsexvSub(prob->Cbar, cell->candidU, rhsx);

	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
		/* only a fraction of subproblems are solved in any iterations */
		if ( scalit(0, 1, &config.RUN_SEED) <= config.SAMPLE_FRACTION ) {
			/* 2. update the right-hand side with exogenous information */
			if ( computeRHS(prob->num, prob->coord, rhsx, cell->omega->vals[obs], cell->candidU, cell->subprob->lp) )
				errMsg("algorithm", "solveSubprobs", "failed to setup subproblem for solve", 0);

			/* 3. solve the subproblem */
			if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
				if ( status == STAT_INFEASIBLE ) {
					/* TODO: subproblem is infeasible add a feasibility cut */
					printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
					writeProblem(cell->subprob->lp, "infeasSub.lp");
					return 1;
				}
				else {
					errMsg("algorithm", "algo", "failed to solve subproblem in solver", 0);
					return 1;
				}
			}

			/* 4. stochastic updates */
			cell->omega->istar[obs] = stocUpdate(cell->subprob->lp, prob->num, prob->coord, prob->Cbar, prob->bBar, cell->sigma, cell->delta,
					cell->omega, rhsx);

#ifdef ALGO_TRACE
			idxSigma = cell->omega->istar[obs];
			idxLambda = cell->sigma->lambdaIdx[idxSigma];
			val = getObjective(cell->subprob->lp, PROB_LP); cummVal += cell->omega->probs[obs]*val;
			printf("\nObs (%d) :: Objective function value = %.3lf\t", obs, val);
			val = cell->sigma->vals[idxSigma].pib + cell->delta->vals[idxLambda][obs].pib -
					vXv(cell->sigma->vals[idxSigma].piC, cell->candidU, prob->coord->colsC, prob->num->cntCcols) -
					vXv(cell->delta->vals[idxLambda][obs].piC, cell->candidU, prob->coord->rvRows+prob->num->rvbOmCnt, prob->num->rvCOmCnt);
			printf("Objective function estimate = %.3lf", val);
#endif
		}
		else
			cell->omega->istar[obs] = -1; /*indicates that the observation was not sampled */
	}

#ifdef ALGO_TRACE
	printf("\n ====> Estimate of expected recourse at Iteration-%d solution = %lf.\n", cell->k, cummVal);
#endif

	mem_free(rhsx);
	return 0;
}//END solveSubprobs()

int computeRHS(numType *num, coordType *coord, vector rhsx, vector observ, vector xk, LPptr lp) {
	sparseVector bomega;
	sparseMatrix Comega;
	vector 	rhs;
	intvec	indices;
	int		cnt, stat1;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val = observ;

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
	Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = observ + num->rvbOmCnt;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)) )
		errMsg("allocation", "solveSubporb", "indices", 0);
	if ( !(rhs = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "solveSubporb", "rhs", 0);

	/* copy right-hand side modified with mean information */
	for ( cnt = 1; cnt <= num->rows; cnt++ ) {
		rhs[cnt] = rhsx[cnt];
		indices[cnt-1] = cnt-1;
	}

	/* change right-hand side with randomness in b */
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* change right-hand side with randomness in transfer matrix */
	rhs = MSparsexvSub(&Comega, xk, rhs);

	/* change the right-hand side in the solver */
	stat1 = changeRHS(lp, num->rows, indices, rhs+1);
	if ( stat1 ) {
		errMsg("solver", "solveSubprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

#ifdef MODIFY_CHECK
	writeProblem(lp, "subprob_k.lp");
#endif
	mem_free(rhs);mem_free(indices);
	return 0;
}//END computeRHS()

int stocUpdate(LPptr lp, numType *num, coordType *coord, sparseMatrix *Cbar, sparseVector *bBar, sigmaType *sigma, deltaType *delta, omegaType *omega,
		vector rhsx) {
	vector	pi;
	double	mubBar;
	int 	idxLambda, idxSigma;
	BOOL	newLambdaFlag;

	/* obtain the dual solution */
	if ( !(pi = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "stocUpdate", "pi", 0);

	if (getDual(lp, pi, num->rows) ) {
		errMsg("solver", "backwardPass", "failed to obtain optimal dual solutions to TDA problem", 0);
		return 1;
	}

	/* compute \bar{\mu} */
	if (calcMu(lp, num->cols, &mubBar) ) {
		errMsg("algorithm", "backwardPass", "failed to compute mu for stochastic updates", 0);
		return 1;
	}

	/* need to calculate new row in the delta matrix if a new pi has been encountered */
	idxLambda = calcDeltaRow(num, coord, omega, delta, pi, &newLambdaFlag);

	/* update the dual information with respect to bBar and Cbar */
	idxSigma = calcSigma(num, coord, bBar, Cbar, pi, mubBar, idxLambda, newLambdaFlag, sigma);

	mem_free(pi);
	return idxSigma;
}//END stocUpdate()

int calcMu(LPptr lp, int numCols, double *mubBar) {
	vector	dj, u;
	intvec	cstat;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);

	if ( getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	/* extra column for eta if the stage problem is a QP */
	if ( !(cstat = (intvec) arr_alloc(numCols+2, int)) )
		errMsg("allocation", "computeMu", "column status", 0);
	if (getBasis(lp, cstat+1, NULL)) {
		errMsg("solver", "computeMu", "failed to get column status", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(cstat); mem_free(dj);

	return 0;
}//END calcMu()

int calcDeltaRow(numType *num, coordType *coord, omegaType *omega, deltaType *delta, vector pi, BOOL *newLambdaFlag) {
	vector 	pixC, newPi;
	int 	obs, n, idxLambda;

	(*newLambdaFlag) = TRUE;

	/* allocate memory to hold duals for rows with random variables */
	if (!(newPi = (vector) arr_alloc(num->rvRowCnt+1, double)) )
		errMsg("allocation", "calcLambda", "current dual vector corresponding to random rows", 0);
	for ( n = 1; n <= num->rvRowCnt; n++ )
		newPi[n] = pi[coord->rvRows[n]];

	/* compare new dual vector with previously observed duals stored in lambda structure */
	idxLambda = 0;
	while (idxLambda < delta->rows) {
		if ( equalVector(newPi, delta->lambda[idxLambda], num->rvRowCnt, config.TOLERANCE) ) {
			(*newLambdaFlag) = FALSE;
			mem_free(newPi);
			return idxLambda;
		}
		idxLambda++;
	}

	delta->rows++;
	delta->lambda[idxLambda] = newPi;

	if ( !(pixC = (vector) arr_alloc(num->rvColCnt, double)) )
		errMsg("allocation", "calcDeltaRow", "pixC", 0);

	/* allocate memory to a new row in delta structure */
	if ( !(delta->vals[idxLambda] = (pixbCType *) arr_alloc(omega->cnt, pixbCType)) )
		errMsg("allocation", "calcDeltaCol", "delta->vals[lambdaIdx]", 0);

	/* For all observations, calculate \pi \times b and \pi \times C */
	for (obs = 0; obs < omega->cnt; obs++) {
		delta->vals[idxLambda][obs].pib = 0;
		for ( n = 1; n <= num->rvRowCnt; n++ )
			delta->vals[idxLambda][obs].pib += delta->lambda[idxLambda][n]*omega->vals[obs][n];
		for ( n = num->rvRowCnt+1; n <= num->numRV; n++ )
			pixC[n-num->rvRowCnt] = delta->lambda[idxLambda][n]*omega->vals[obs][n];
		delta->vals[idxLambda][obs].piC = pixC;
	}

	return idxLambda;
}//END calcDeltaRow()

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar, int idxLambda, BOOL newLambdaFlag,
		sigmaType *sigma) {
	double pibBar;
	vector	piCBar, temp;
	int		cnt;

	pibBar = vXvSparse(pi, bBar) + mubBar;

	temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->colsC, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= config.TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, config.TOLERANCE)) {
					mem_free(piCBar);
					return cnt;
				}
			}
		}
	}

	sigma->vals[sigma->cnt].pib  = pibBar;
	sigma->vals[sigma->cnt].piC  = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;

	return sigma->cnt++;
}//END calcSigma()
