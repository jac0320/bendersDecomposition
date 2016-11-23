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


/* benders.c */
void parseCmdLine(int argc, string *argv, string probName);
int readConfig(int argc, string probName);


#endif /* BENDERS_H_ */
