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
	double objVal, objEst;

	if ( cell->k <= config.MIN_ITER)
		return FALSE;

	if ( !(primal = (vector) arr_alloc(cell->master->macsz+1, double)) )
		errMsg("allocation", "optimal", "primal", 0);
	objVal = vXv(cell->master->objx-1, cell->candidU, NULL, cell->master->macsz);
	if (cell->master->type == PROB_LP)
		objEst = cell->candidEst;
	else if (cell->master->type == PROB_QP)
		objEst = cell->incumbEst;

	if ( DBL_ABS(objVal - objEst) <= config.TOLERANCE ) {
		printf("\nOptimal objective function value = %lf\n", objVal);
		return TRUE;
	}
	else
		printf("Obective value = %lf; Incumbent estimate = %lf\n", objVal, objEst);
//		printf("+"); fflush(stdout);

	return FALSE;
}//END optimal()
