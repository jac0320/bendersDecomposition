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

#undef INPUT_CHECK
#undef MODIFY_CHECK
#define ALGO_TRACE

typedef struct{
	int			SAA;
	int 		SAA_OBS;
	long long	SAA_SEED;
	int 		MAX_ITER;
	int			MULTICUT;
	int 		PROXIMAL;
	int			MASTER_TYPE;
	long long	RUN_SEED;
	double		SAMPLE_FRACTION;
	double		TOLERANCE;
}configType;

/* Omega stores the set of observations. Each observation consists of a vector of realization of random variables with discrete distributions and
 * the associated probability. */
typedef struct {
	int 	cnt;
	vector	probs;
	vector	*vals;
	int     numRV;
	intvec	istar;
} omegaType;

typedef struct {
	double	pib;
	vector  piC;
}pixbCType;

typedef struct{
	int 		cnt;
	pixbCType	*vals;
	intvec		lambdaIdx;
}sigmaType;

typedef struct {
	int			rows;
	int			cols;
	vector  	*lambda;
	pixbCType	**vals;
}deltaType;

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
	int			k;
	int			maxCuts;
	oneProblem 	*master;
	oneProblem 	*subprob;
	omegaType	*omega;
	sigmaType	*sigma;
	deltaType	*delta;
	cutsType	*cuts;
	vector		candidU;
	vector		incumbU;
}cellType;

/* benders.c */
void parseCmdLine(int argc, string *argv, string probName);
int readConfig(int argc, string probName);

/* algo.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(probType **prob, stocType *stoc, vector u0);
oneProblem *newMaster(probType *prob, cutsType *cuts, omegaType *omega);
oneProblem *newSubproblem(probType *prob);
omegaType *newOmega(stocType *stoc);
sigmaType *newSigma(int numIter);
deltaType *newDelta(int numIter, int maxObs);
void freeCellType(cellType *cell);
void freeOmegaType(omegaType *omega);
void freeSigmaType(sigmaType *sigma);
void freeDeltaType (deltaType *delta);

/* subproblem.c */
int solveSubprobs(probType *prob, cellType *cell);
int computeRHS(numType *num, coordType *coord, vector rhsx, vector observ, vector U, LPptr lp);
int stocUpdate(LPptr lp, numType *num, coordType *coord, sparseMatrix *Cbar, sparseVector *bBar, sigmaType *sigma, deltaType *delta, omegaType *omega,
		vector rhsx);
int calcMu(LPptr lp, int numCols, double *mubBar);
int calcDeltaRow(numType *num, coordType *coord, omegaType *omega, deltaType *delta, vector pi, BOOL *newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar, int idxLambda, BOOL newLambdaFlag,
		sigmaType *sigma);

/* cuts.c */
int formSingleCut(probType **prob, cellType *cell);
int formMultiCut(probType **prob, cellType *cell);
oneCut *newCut(int betaLen);
int addCut(LPptr lp, cutsType *cuts, int numRows, int betaLen, intvec betaIdx, int maxCuts, int etaIndex, oneCut *cut);
void freeCutsType(cutsType *cuts);
void freeOneCut(oneCut *cut);

#endif /* BENDERS_H_ */
