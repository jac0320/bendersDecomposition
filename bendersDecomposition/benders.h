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
	int			SAA;
	int 		SAA_OBS;
	long long	SAA_SEED;
	int 		MAX_ITER;
	int			MULTICUT;
	int 		PROXIMAL;
	int			MASTER_TYPE;
	long long	RUN_SEED;
}configType;

/* Omega stores the set of observations. Each observation consists of a vector of realization of random variables with discrete distributions and
 * the associated probability. */
typedef struct {
	int 	cnt;
	vector	probs;
	vector	*vals;
	int     numRV;
} omegaType;


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
	omegaType	*omega;
	cutsType	*cuts;
}cellType;

/* benders.c */
void parseCmdLine(int argc, string *argv, string probName);
int readConfig(int argc, string probName);

/* algo.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(probType **prob, stocType *stoc);
oneProblem *newMaster(probType *prob, cutsType *cuts, omegaType *omega);
oneProblem *newSubproblem(probType *prob);
omegaType *newOmega(stocType *stoc);
void freeCellType(cellType *cell);
void freeOmegaType(omegaType *omega);

#endif /* BENDERS_H_ */
