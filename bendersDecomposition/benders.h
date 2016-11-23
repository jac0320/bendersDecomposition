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

typedef struct{
	int MAX_ITER;
}configType;

typedef struct {

}cellType;

/* benders.c */
void parseCmdLine(int argc, string *argv, string probName);
int readConfig(int argc, string probName);

/* setup.c */
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(probType **prob);
void freeCellType(cellType *cell);

#endif /* BENDERS_H_ */
