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

	cell->master = newMaster(prob[0], NULL);
	cell->subprob = newSubproblem(prob[1]);

	return cell;
}//END newCell()

oneProblem *newMaster(probType *prob, cutsType *cuts) {
	oneProblem *master = NULL;
	char		*q, cutName[NAMESIZE] = { "Old    " };
	int			j, cnt, idx, rowOffset, colOffset, cutsCnt = 0;

	if ( cuts != NULL )
		cutsCnt = cuts->cnt;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("allocation", "newMaster", "copy", 0);
	master->bdl  = NULL; master->bdu 	 = NULL; master->cname  = NULL; master->cstore = NULL; master->ctype  = NULL;
	master->lp   = NULL; master->matbeg  = NULL; master->matcnt = NULL; master->matind = NULL; master->matval = NULL;
	master->name = NULL; master->objname = NULL; master->objx   = NULL; master->rhsx   = NULL; master->rname  = NULL;
	master->rstore=NULL; master->senx    = NULL;

	/* assign values to scalar terms */
	master->objsen = prob->sp->objsen;
	master->marsz = prob->sp->mar + cutsCnt;	master->mar = prob->sp->mar + cutsCnt;	master->mac = prob->sp->mac; /* Note: mac does not include eta columms */
	if ( config.MULTICUT ) {
		// TODO: figure out how many columns to add
		/* S additional columns to hold eta's (surrogates for individual recourse functions */
		master->macsz = prob->sp->mac; master->cstorsz = prob->sp->cstorsz + NAMESIZE;
	}
	else {
		/* one additional column to hold eta (surrogate for expected recourse function */
		master->macsz = prob->sp->mac+1; master->cstorsz = prob->sp->cstorsz + NAMESIZE;
	}
	master->matsz  = prob->sp->matsz + (cutsCnt)*(master->mac + 1);
	master->numnz  = prob->sp->numnz;
	master->rstorsz = prob->sp->rstorsz + NAMESIZE*(cutsCnt);
	master->numInt = prob->sp->numInt; master->numBin = prob->sp->numBin;
	if ( master->numInt + master->numBin > 0 ) {
		if ( config.PROXIMAL == 0) {
			if ( config.MASTER_TYPE == PROB_LP || config.MASTER_TYPE == PROB_QP)
				printf("Warning :: Ignoring integrality and solving a relaxed problem.\n");
			master->type = config.MASTER_TYPE;
		}
		else {
			if ( config.MASTER_TYPE == PROB_MILP || config.MASTER_TYPE == PROB_MIQP ) {
				errMsg("configuration", "newMaster", "No regularization supported for MILPs and MIQPs", 0);
				return NULL;
			}
			printf("Warning :: Ignoring integrality and solving a regularized, relaxed problem.\n");
			master->type = PROB_QP;
		}
	}
	else
		(config.PROXIMAL == 1) ? (master->type = PROB_QP): (master->type = config.MASTER_TYPE);

	/* Make all allocations of known sizes, as calculated above */
	if (!(master->name = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newMaster", "master->name",0);
	if (!(master->objname = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newMaster", "master->objname",0);
	if (!(master->objx = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->objx",0);
	if (!(master->bdl = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->bdl",0);
	if (!(master->bdu = arr_alloc(master->macsz, double)))
		errMsg("allocation", "newMaster", "master->bdu",0);
	if (!(master->rhsx = arr_alloc(master->marsz, double)))
		errMsg("allocation", "newMaster", "master->rhsx",0);
	if (!(master->senx = arr_alloc(master->marsz, char)))
		errMsg("allocation", "newMaster", "master->senx",0);
	if (!(master->matbeg = arr_alloc(master->macsz, int)))
		errMsg("allocation", "newMaster", "master->matbeg",0);
	if (!(master->matcnt = arr_alloc(master->macsz, int)))
		errMsg("allocation", "newMaster", "master->matcnt",0);
	if (!(master->cname = arr_alloc(master->macsz, string)))
		errMsg("allocation", "newMaster", "master->cname",0);
	if (!(master->cstore = arr_alloc(master->cstorsz, char)))
		errMsg("allocation", "newMaster", "master->cstore",0);
	if (!(master->rname = arr_alloc(master->marsz, string)))
		errMsg("allocation", "newMaster", "master->rname",0);
	if (!(master->rstore = arr_alloc(master->rstorsz, char)))
		errMsg("allocation", "newMaster", "master->rstore",0);
	if (!(master->matval = arr_alloc(master->matsz, double)))
		errMsg("allocation", "newMaster", "master->matval",0);
	if (!(master->matind = arr_alloc(master->matsz, int)))
		errMsg("allocation", "newMaster", "master->matind",0);
	if (!(master->ctype = arr_alloc(master->macsz, char)))
		errMsg("allocation", "newMaster", "master->ctype",0);

	/* First copy information directly from the original subproblem. */
	/* Copy the subproblem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	j = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		master->cstore[j++] = *q;

	j = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		master->rstore[j++] = *q;

	strcpy(master->name, prob->sp->name);
	strcpy(master->objname, prob->sp->objname);

	/* Calculate difference in pointers for oneProblems in probType and cellType for row and column names */
	colOffset = master->cstore - prob->sp->cname[0];
	rowOffset = master->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original subproblem */
	cnt = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		master->objx[j] = prob->sp->objx[j];
		master->ctype[j] = prob->sp->ctype[j];
		master->bdu[j] = prob->sp->bdu[j];
		master->bdl[j] = prob->sp->bdl[j];
		master->cname[j] = prob->sp->cname[j] + colOffset;
		master->matbeg[j] = cnt;
		master->matcnt[j] = prob->sp->matcnt[j];
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			master->matval[cnt] = prob->sp->matval[idx];
			master->matind[cnt] = prob->sp->matind[idx];
			cnt++;
		}
	}

	/* Copy all row information from the original subproblem */
	for (j = 0; j < prob->sp->mar; j++) {
		master->rhsx[j] = prob->sp->rhsx[j];
		master->senx[j] = prob->sp->senx[j];
		master->rname[j] = prob->sp->rname[j] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	if (config.MULTICUT == 0) {
		strcpy(master->cstore + master->cstorsz, "eta");
		master->cname[master->mac] = master->cstore + master->cstorsz;
		master->objx[master->mac] = 1.0;
		master->ctype[master->mac] = 'C';
		master->bdu[master->mac] = INFBOUND;
		master->bdl[master->mac] = prob->lb;
		master->matbeg[master->mac] = cnt;
		master->matcnt[master->mac] = cutsCnt;
	}
	else {
		printf("Multicut version is not supported\n");
		return NULL;
	}

	if (cutsCnt > 0) {
		/* Now copy information from the cuts into the new master problem. */
		for (idx = 0; idx < cuts->cnt; idx++) {
			master->matval[idx + cnt] = 1.0;
			master->matind[idx + cnt] = master->mar + idx;
		}
		/* Copy the constraint coefficients */
		for (j = 0; j < master->mac; j++) {
			cnt = master->matbeg[j] + master->matcnt[j];
			master->matcnt[j] += cuts->cnt;

			for (idx = 0; idx < cuts->cnt; idx++) {
				master->matval[cnt + idx] = cuts->vals[idx]->beta[j + 1];
				master->matind[cnt + idx] = master->mar + idx;
			}
		}

		/* Give names to the cut constraints, and add rhs values */
		idx = master->rstorsz;
		for (cnt = 0; cnt < cuts->cnt; cnt++) {
			cutName[3] = '0' + cnt / 1000 % 10;
			cutName[4] = '0' + cnt / 100 % 10;
			cutName[5] = '0' + cnt / 10 % 10;
			cutName[6] = '0' + cnt / 1 % 10;
			master->rname[master->mar + cnt] = master->rstore + idx;
			strcpy(master->rname[master->mar + cnt], cutName);
			idx += strlen(cutName) + 1;
			master->rhsx[master->mar + cnt] = cuts->vals[cnt]->alpha;
			master->senx[master->mar + cnt] = 'G';
		}
	}

	/* Load the copy into CPLEX now */
	master->lp = setupProblem(master->name, master->type, master->macsz, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,
			master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "newMaster", "failed to setup master problem in the solver",0);
		return NULL;
	}

#ifdef INPUT_CHECK
	writeProblem(master->lp, "master.lp");
#endif

	return master;
}//END newMaster()

oneProblem *newSubproblem(probType *prob) {
	oneProblem  *subprob = NULL;
	char		*q;
	int			j, cnt, idx, rowOffset, colOffset;

	if (!(subprob = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("allocation", "newSubproblem", "copy", 0);
	subprob->bdl  = NULL; subprob->bdu 	  = NULL; subprob->cname  = NULL; subprob->cstore = NULL; subprob->ctype = NULL;
	subprob->lp   = NULL; subprob->matbeg = NULL; subprob->matcnt = NULL; subprob->matind = NULL; subprob->matval = NULL;
	subprob->name = NULL; subprob->objname= NULL; subprob->objx   = NULL; subprob->rhsx	  = NULL; subprob->rname = NULL;
	subprob->rstore=NULL; subprob->senx = NULL;

	/* assign values to scalar terms */
	subprob->objsen = prob->sp->objsen;
	subprob->matsz = prob->sp->matsz;		subprob->numnz = prob->sp->numnz;
	subprob->marsz = prob->sp->mar;			subprob->mar = prob->sp->mar;
	subprob->macsz = prob->sp->mac; 		subprob->mac = prob->sp->mac;
	subprob->cstorsz = prob->sp->cstorsz;	subprob->rstorsz = prob->sp->rstorsz;
	subprob->numInt = prob->sp->numInt;		subprob->numBin = prob->sp->numBin;
	if ( subprob->numInt + subprob->numBin > 0 ) {
		subprob->type = PROB_MILP;
		errMsg("read", "newSubproblem", "subproblem has discrete variables", 0);
		return NULL;
	}
	else
		subprob->type = PROB_LP;

	/* Make all allocations of known sizes, as calculated above */
	if (!(subprob->name = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newSubproblem", "subprob->name",0);
	if (!(subprob->objname = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "newSubproblem", "subprob->objname",0);
	if (!(subprob->objx = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->objx",0);
	if (!(subprob->bdl = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->bdl",0);
	if (!(subprob->bdu = arr_alloc(subprob->macsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->bdu",0);
	if (!(subprob->rhsx = arr_alloc(subprob->marsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->rhsx",0);
	if (!(subprob->senx = arr_alloc(subprob->marsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->senx",0);
	if (!(subprob->matbeg = arr_alloc(subprob->macsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matbeg",0);
	if (!(subprob->matcnt = arr_alloc(subprob->macsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matcnt",0);
	if (!(subprob->cname = arr_alloc(subprob->macsz, string)))
		errMsg("allocation", "newSubproblem", "subprob->cname",0);
	if (!(subprob->cstore = arr_alloc(subprob->cstorsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->cstore",0);
	if (!(subprob->rname = arr_alloc(subprob->marsz, string)))
		errMsg("allocation", "newSubproblem", "subprob->rname",0);
	if (!(subprob->rstore = arr_alloc(subprob->rstorsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->rstore",0);
	if (!(subprob->matval = arr_alloc(subprob->matsz, double)))
		errMsg("allocation", "newSubproblem", "subprob->matval",0);
	if (!(subprob->matind = arr_alloc(subprob->matsz, int)))
		errMsg("allocation", "newSubproblem", "subprob->matind",0);
	if (!(subprob->ctype = arr_alloc(subprob->macsz, char)))
		errMsg("allocation", "newSubproblem", "subprob->ctype",0);

	/* First copy information directly from the original subproblem. */
	/* Copy the subproblem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	j = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		subprob->cstore[j++] = *q;

	j = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		subprob->rstore[j++] = *q;

	strcpy(subprob->name, prob->sp->name);
	strcpy(subprob->objname, prob->sp->objname);

	/* Calculate difference in pointers for oneProblems in probType and cellType for row and column names */
	colOffset = subprob->cstore - prob->sp->cname[0];
	rowOffset = subprob->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original subproblem */
	cnt = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		subprob->objx[j] = prob->sp->objx[j];
		subprob->ctype[j] = prob->sp->ctype[j];
		subprob->bdu[j] = prob->sp->bdu[j];
		subprob->bdl[j] = prob->sp->bdl[j];
		subprob->cname[j] = prob->sp->cname[j] + colOffset;
		subprob->matbeg[j] = cnt;
		subprob->matcnt[j] = prob->sp->matcnt[j];
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			subprob->matval[cnt] = prob->sp->matval[idx];
			subprob->matind[cnt] = prob->sp->matind[idx];
			cnt++;
		}
	}

	/* Copy all row information from the original subproblem */
	for (j = 0; j < prob->sp->mar; j++) {
		subprob->rhsx[j] = prob->sp->rhsx[j];
		subprob->senx[j] = prob->sp->senx[j];
		subprob->rname[j] = prob->sp->rname[j] + rowOffset;
	}

	/* load the information on to the solver */
	subprob->lp = setupProblem(subprob->name, subprob->type, subprob->mac, subprob->mar, subprob->objsen, subprob->objx, subprob->rhsx, subprob->senx,
			subprob->matbeg, subprob->matcnt, subprob->matind, subprob->matval, subprob->bdl, subprob->bdu, NULL, subprob->cname, subprob->rname,
			subprob->ctype);
	if ( subprob->lp == NULL ) {
		errMsg("solver", "newSubprob", "subprob",0);
		return NULL;
	}

#ifdef INPUT_CHECK
	writeProblem(subprob->lp, "subproblem.lp");
#endif

	return subprob;
}//END newSubproblem()

void freeCellType(cellType *cell) {

	if (cell) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->subprob) freeOneProblem(cell->subprob);
		mem_free(cell);
	}

}//END freeCellType()
