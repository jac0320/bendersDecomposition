/*
 * benders.h
 *
 *  Created on: Nov 22, 2016
 *      Author: Harsha Gangammanavar
 *  Instituion: Southern Methodist University 
 *      e-mail: harsha(at)smu(dot)edu 
 */

#ifndef BENDERS_H_
#define BENDERS_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#define INPUT_CHECK

typedef struct{
	int MAX_ITER;
	int	MULTICUT;
	int PROXIMAL;
	int	MASTER_TYPE;
}configType;

typedef struct {
	double 	alpha;					/* scalar value for the right-hand side */
	vector 	beta;					/* coefficients for the master program's primal variables */
	int 	rowNum;					/* row number for master problem in solver */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	char  	sense;					/* sense of the cut being added */
}oneCut;

typedef struct{
	int 	cnt;
	oneCut	**vals;
}cutsType;

typedef struct {
	oneProblem 	*master;
	oneProblem 	*subprob;
	cutsType	*cuts;
}cellType;

/* benders.c */
void parseCmdLine(int argc, string *argv, string probName);
int readConfig(int argc, string probName);

/* algo.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(probType **prob);
oneProblem *newMaster(probType *prob, cutsType *cuts);
oneProblem *newSubproblem(probType *prob);
void freeCellType(cellType *cell);

#endif /* BENDERS_H_ */
