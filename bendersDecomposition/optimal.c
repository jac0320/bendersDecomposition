/*
 * optimal.c
 *
 *  Created on: Nov 29, 2016
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University.
 *  	Report: harsha(at)smu(dot)edu
 */

#include "benders.h"

configType config;

BOOL optimal(probType *prob, cellType *cell) {
	vector primal;
	double objVal;
	int    n;

	if ( cell->k <= config.MIN_ITER)
		return FALSE;

	/* obtain primal solution */
	if ( !(primal = (vector) arr_alloc(cell->master->macsz+1, double)) )
		errMsg("allocation", "optimal", "primal", 0);
	if ( getPrimal(cell->master->lp, primal, cell->master->macsz) ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the optimal primal solution", 0);
		return 1;
	}
	objVal = vXv(cell->master->objx-1, cell->candidU, NULL, prob->num->cols);
	for ( n = prob->num->cols; n < cell->master->macsz; n++ )
		objVal += primal[n+1];

	if (cell->master->type == PROB_LP) {
		if ( DBL_ABS(objVal - cell->candidEst) <= config.TOLERANCE ) {
			printf("\nOptimal objective function value = %lf\n", objVal);
			return TRUE;
		}
	}
	else if (cell->master->type == PROB_QP) {
		if ( DBL_ABS(objVal - cell->incumbEst) <= config.TOLERANCE ) {
			printf("\nOptimal objective function value = %lf\n", objVal);
			return TRUE;
		}
	}

	return FALSE;
}//END optimal()
