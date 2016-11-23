/*
 * algo.c
 *
 *  Created on: Nov 22, 2016
 *      Author: Harsha Gangammanavar
 *  Instituion: Southern Methodist University 
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

	TERMINATE:
	freeCellType(cell);
	freeProbType(prob, tim->numRows);
	return 0;
}//END algo()

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell) {
	vector	xk = NULL, lb = NULL;

	/* setup mean value problem which will act as reference for all future computations */
	xk = meanProblem(orig, stoc);
	if ( xk == NULL ) {
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
	(*cell) = newCell((*prob));
	if ((*cell) == NULL) {
		errMsg("setup", "setupAlgo", "failed to setup the cell used by the algorithm", 0);
		goto TERMINATE;
	}

	mem_free(lb);
	mem_free(xk);
	return 0;

	TERMINATE:
	if (lb) mem_free(lb);
	if (xk) mem_free(xk);
	return 1;
}//END setupAlgo()

cellType *newCell(probType **prob) {
	cellType *cell;

	if (!(cell = (cellType *) mem_malloc (sizeof(cellType))))
		errMsg("allocation", "newCell", "cell",0);

	return cell;
}//END newCell()

void freeCellType(cellType *cell) {

	if (cell) {
		mem_free(cell);
	}

}//END freeCellType()
