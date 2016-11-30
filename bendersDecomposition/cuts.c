/*
 * cuts.c
 *
 *  Created on: Nov 24, 2016
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *      e-mail: harsha(at)smu(dot)edu 
 */

#include "benders.h"

configType config;

int formSingleCut(probType **prob, cellType *cell) {
	oneCut	*cut;
	vector 	pixC;
	double	arg, argmax;
	int 	obs, cnt, i, deltaRow, istar;

	/* 0. allocate memory to new cut */
	cut = newCut(prob[1]->num->cntCcols);
	if ( cut == NULL ) {
		errMsg("algorithm", "formSingleCut", "failed to allocate memory to the new cut", 0);
		return -1;
	}

	if ( !(pixC = (vector) arr_alloc(cell->sigma->cnt, double)) )
		errMsg("allocation", "formSingleCut", "pixC", 0);

	/* calculate [\bar{b}^\top \pi - (\bar{C}_t^\top \pi)x] one at a time for all entries in sigma */
	for ( cnt = 0; cnt < cell->sigma->cnt; cnt++ ) {
		for (i = 1; i <= prob[1]->num->cntCcols; i++)
			pixC[cnt] += cell->sigma->vals[cnt].piC[i] * cell->candidU[prob[1]->coord->colsC[i]];
		pixC[cnt] = cell->sigma->vals[cnt].pib - pixC[cnt];
	}

	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
		if ( cell->omega->istar[obs] < 0 ) {
			/* 1. Identify the best dual solution for observations which was not selected during sampling */
			argmax = -DBL_MAX;
			for ( cnt = 0; cnt < cell->sigma->cnt; cnt++ ) {
				deltaRow = cell->sigma->lambdaIdx[cnt];

				/* Start with (\pi^\top \tilde{\Delta b}) + [\bar{b}^\top \pi - (\bar{C}_t^\top \pi)x] */
				arg = pixC[cnt] + cell->delta->vals[deltaRow][obs].pib;
				/* Subtract (\tilde{C}_t^\top \pi)*u_t */
				for (i = 1; i <= prob[1]->num->rvColCnt; i++)
					arg -= cell->delta->vals[deltaRow][cnt].piC[i] * cell->candidU[prob[1]->coord->rvCols[i]];

				if (arg > argmax) {
					argmax = arg;
					cell->omega->istar[obs] = cnt;
				}
			}
		}

		/* 2a. Average using these pi's to calculate the cut coefficients. Here we are multiplying the coefficients by the probability
			   of each observations: Intercept term first */
		istar = cell->omega->istar[obs];
		cut->alpha += (cell->sigma->vals[istar].pib + cell->delta->vals[cell->sigma->lambdaIdx[istar]][obs].pib)* cell->omega->probs[obs];

		/* 2b. Slope term next */
		for (i = 1; i <= prob[1]->num->cntCcols; i++)
			cut->beta[prob[1]->coord->colsC[i]] += cell->sigma->vals[istar].piC[i] * cell->omega->probs[obs];
		for (i = 1; i <= prob[1]->num->rvColCnt; i++)
			cut->beta[prob[1]->coord->rvCols[i]] += cell->delta->vals[cell->sigma->lambdaIdx[istar]][obs].piC[i] * cell->omega->probs[obs];
	}

	mem_free(pixC);

	/* 3. Check to see if the cut is already present in optimality cuts pool */
	cnt = 0;
	if ( cnt < cell->cuts->cnt ) {
		if ( DBL_ABS(cell->cuts->vals[cnt]->alpha - cut->alpha) <= config.TOLERANCE )
			if (equalVector(cell->cuts->vals[cnt]->beta, cut->beta, prob[1]->num->cntCcols, config.TOLERANCE)) {
				/* the cut already exists in the pool */
				freeOneCut(cut);
				return 0;
			}
		cnt++;
	}

	if ( config.PROXIMAL == 1 )
		cut->alphaIncumb = cut->alpha + vXv(cut->beta, cell->incumbU, prob[1]->coord->colsC, prob[1]->num->cntCcols);

	/* 4. Add cut to the master problem cut structure as well as on the solver */
	istar = addCut(cell->master->lp, cell->cuts, prob[0]->num->rows, prob[1]->num->cntCcols, prob[1]->coord->colsC,
					cell->maxCuts, prob[0]->num->cols+0, cut);
	if ( istar < 0 ) {
		errMsg("algorithm", "formSingleCut", "failed to add the cut to master problem", 0);
		freeOneCut(cut);
		return -1;
	}

	cell->candidEst = vXv(cell->master->objx-1, cell->candidU, NULL, prob[0]->num->cols) +
			cut->alpha - vXv(cut->beta, cell->candidU, prob[1]->coord->colsC, prob[1]->num->cntCcols);

#ifdef ALGO_TRACE
	printf(" ====> Cut height at Iteration-%d solution                    = %lf\n", cell->k,
			cell->candidEst);
#endif

	return 0;
}//END formSingleCut()

int formMultiCut(probType **prob, cellType *cell) {
	oneCut	*cut;
	vector 	pixC;
	double	arg, argmax;
	int 	obs, cnt, i, deltaRow, istar;


	if ( !(pixC = (vector) arr_alloc(cell->sigma->cnt, double)) )
		errMsg("allocation", "formMultiCut", "pixC", 0);

	/* calculate [\bar{b}^\top \pi - (\bar{C}_t^\top \pi)x] one at a time for all entries in sigma */
	for ( cnt = 0; cnt < cell->sigma->cnt; cnt++ ) {
		for (i = 1; i <= prob[1]->num->cntCcols; i++)
			pixC[cnt] += cell->sigma->vals[cnt].piC[i] * cell->candidU[prob[1]->coord->colsC[i]];
		pixC[cnt] = cell->sigma->vals[cnt].pib - pixC[cnt];
	}

	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
		if ( cell->omega->istar[obs] < 0 ) {
			/* 1. Identify the best dual solution for observations which was not selected during sampling */
			argmax = -DBL_MAX;
			for ( cnt = 0; cnt < cell->sigma->cnt; cnt++ ) {
				deltaRow = cell->sigma->lambdaIdx[cnt];

				/* Start with (\pi^\top \tilde{\Delta b}) + [\bar{b}^\top \pi - (\bar{C}_t^\top \pi)x] */
				arg = pixC[cnt] + cell->delta->vals[deltaRow][obs].pib;
				/* Subtract (\tilde{C}_t^\top \pi)*u_t */
				for (i = 1; i <= prob[1]->num->rvColCnt; i++)
					arg -= cell->delta->vals[deltaRow][cnt].piC[i] * cell->candidU[prob[1]->coord->rvCols[i]];

				if (arg > argmax) {
					argmax = arg;
					cell->omega->istar[obs] = cnt;
				}
			}
		}

		/* 2. allocate memory to new cut */
		cut = newCut(prob[1]->num->cntCcols);
		if ( cut == NULL ) {
			errMsg("algorithm", "formMultiCut", "failed to allocate memory to the new cut", 0);
			return -1;
		}

		/* 3a. Average using these pi's to calculate the cut coefficients. Here we are multiplying the coefficients by the probability
			   of each observations: Intercept term first */
		istar = cell->omega->istar[obs];
		cut->alpha += cell->sigma->vals[istar].pib + cell->delta->vals[istar][obs].pib;

		/* 3b. Slope term next */
		for (i = 1; i <= prob[1]->num->cntCcols; i++)
			cut->beta[prob[1]->coord->colsC[i]] += cell->sigma->vals[istar].piC[i];
		for (i = 1; i <= prob[1]->num->rvColCnt; i++)
			cut->beta[prob[1]->coord->rvCols[i]] += cell->delta->vals[istar][obs].piC[i];

		/* 4. Check to see if the cut is already present in optimality cuts pool */
		cnt = 0;
		if ( cnt < cell->cuts->cnt ) {
			if ( cell->cuts->vals[cnt]->beta[0] == prob[1]->num->cols+obs)
				if ( DBL_ABS(cell->cuts->vals[cnt]->alpha - cut->alpha) <= config.TOLERANCE )
					if (equalVector(cell->cuts->vals[cnt]->beta, cut->beta, prob[0]->num->cols, config.TOLERANCE)) {
						/* the cut already exists in the pool */
						freeOneCut(cut);
						return 0;
					}
			cnt++;
		}
		cut->alphaIncumb = cut->alpha + vXv(cut->beta, cell->incumbU, NULL, prob[0]->num->cols);

		/* 5. Add cut to the master problem cut structure as well as on the solver */
		istar = addCut(cell->master->lp, cell->cuts, prob[0]->num->rows, prob[1]->num->cntCcols, prob[1]->coord->colsC,
				cell->maxCuts, prob[0]->num->cols+obs, cut);
		if ( istar < 0 ) {
			errMsg("algorithm", "formMultiCut", "failed to add the cut stage problem", 0);
			return -1;
		}

#ifdef ALGO_TRACE
	printf(" ====> Cut height at Iteration-%d solution for observation-%d =  %lf\n", cell->k, obs,
			cut->alpha - vXv(cut->beta, cell->candidU, prob[1]->coord->colsC, prob[1]->num->cntCcols));
#endif
	}

	mem_free(pixC);
	return 0;
}//END formMultiCut()

/* subroutine to the allocate memory to oneCut structure and initialize its elements with default values */
oneCut *newCut(int betaLen) {
	oneCut *cut;

	if ( !(cut = (oneCut *) mem_malloc (sizeof(oneCut))) )
		errMsg("allocation", "newCut", "cut", 0);

	if ( !(cut->beta = (vector) arr_alloc(betaLen+1, double)) )
		errMsg("allocation", "newCut", "cut->beta", 0);
	cut->alpha 	 = 0.0;
	cut->beta[0] = 1.0;
	cut->rowNum = -1;
	cut->alphaIncumb = 0.0;
	cut->sense = 'G';

	return cut;
}//END newCut()

int addCut(LPptr lp, cutsType *cuts, int numRows, int betaLen, intvec betaIdx, int maxCuts, int etaIndex, oneCut *cut) {
	intvec	indices;
	int		n;

	if ( !(indices = (intvec) arr_alloc(betaLen+1, int)) )
		errMsg("allocation", "addCut", "indices", 0);
	for ( n = 1; n <= betaLen; n++ )
		indices[n] = betaIdx[n]-1;
	indices[0] = etaIndex;

	/* TODO: make sure there is room to add a new cut */
	if (cuts->cnt >= maxCuts) {
		errMsg("algorithm", "addCut", "ran out of memory for cuts", 0);
		goto TERMINATE;
	}

	if ( addRow(lp, betaLen+1, cut->alpha, cut->sense, 0, indices, cut->beta) ) {
		errMsg("solver", "addCut", "failed to add the new cut to solver problem", 0);
		goto TERMINATE;
	}

	// TODO: verify the cut row number once you switch to quadratics
	cut->rowNum = numRows + cuts->cnt;
	cuts->vals[cuts->cnt] = cut;

	mem_free(indices);
	return cuts->cnt++;

	TERMINATE:
	mem_free(indices);
	return -1;
}//END addCut()

void freeCutsType(cutsType *cuts) {
	int n;

	if (cuts->vals) {
		for ( n = 0; n < cuts->cnt; n++ )
			if (cuts->vals[n]) freeOneCut(cuts->vals[n]);
		mem_free(cuts->vals);
	}
	mem_free(cuts);

}//END freeCutsType()

void freeOneCut(oneCut *cut) {

	if (cut->beta) mem_free(cut->beta);
	mem_free(cut);

}//END freeOneCut()
