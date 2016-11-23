/*
 * benders.c
 *
 *  Created on: Nov 22, 2016
 *      Author: Harsha Gangammanavar
 *  Instituion: Southern Methodist University 
 *      e-mail: harsha(at)smu(dot)edu 
 */

#include "benders.h"

long 		MEM_USED = 0;	/* Amount of memory allocated each iteration */
string   	outputDir;		/* output directory */
configType	config;			/* configuration structure */

int main(int argc, char *argv[]) {
	oneProblem 	*orig;
	timeType	*tim;
	stocType	*stoc;
	char probName[NAMESIZE], probPath[BLOCKSIZE], buffer[2*BLOCKSIZE];

	/* open solver environment */
	openSolver();

	/* request for problem name to be solved, the default path is provided in the configuration file */
	parseCmdLine(argc, argv, probName);

	/* read algorithm configuration file */
	if ( readConfig(argc, probPath) ) {
		errMsg("read", "main", "failed to read configuration file for the algorithm", 0);
		goto TERMINATE;
	}

	/* read problem files */
	if ( readFiles(probPath, probName, &orig, &tim, &stoc) ) {
		errMsg("read", "main", "failed to read problem files", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in config file and the input problem name */
	strcat(outputDir,probName);
	strcat(outputDir,"/");
	sprintf(buffer, "mkdir %s", outputDir);
	system(buffer);
	strcat(outputDir,"benders/");
	sprintf(buffer, "mkdir %s", outputDir);
	if (system(buffer)){
		sprintf(buffer, "rm -r %s*", outputDir);
		system(buffer);
	}

	/* launch the algorithm */
	if (algo(orig, stoc, tim) ) {
		errMsg("algorithm", "main", "failed to solve the problem using integer-SD", 0);
		goto TERMINATE;
	}
	printf("\n--------------------------------------------------------------------------------------------------------------------------------------");
	printf("\nSuccessfully solved the problem (%s) with stochastic Benders decomposition algorithm. Solution files are written to output directory.\n", orig->name);
	printf("--------------------------------------------------------------------------------------------------------------------------------------");

	/* close solver environment and release all structures */
	TERMINATE:
	mem_free(outputDir);
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();

	return 0;
}//END main()

void parseCmdLine(int argc, string *argv, string probName) {

	switch (argc) {
	case 2:
		strcpy(probName, argv[1]);
		break;
	case 1:
		printf("Please enter the name of the problem: ");
		scanf("%s", probName);
		break;
	default:
		break;
	}

}//END parseCmdLine

int readConfig(int argc, string probPath) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status;

	fptr = fopen("./benders.conf", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	outputDir = (string) mem_malloc(BLOCKSIZE*sizeof(char));

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "INPUTDIR")))
			fscanf(fptr, "%s", probPath);
		else if (!(strcmp(line, "OUTPUTDIR")))
			fscanf(fptr, "%s", outputDir);
		else if (!(strcmp(line, "MULTICUT")))
			fscanf(fptr, "%d", &config.MULTICUT);
		else if (!(strcmp(line, "PROXIMAL")))
			fscanf(fptr, "%d", &config.PROXIMAL);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	return 0;
}//END readConfig()
