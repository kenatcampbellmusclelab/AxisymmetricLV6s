// runModel - version for parameter estimation
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>

void initialgeometry(int Nmu, int Nnu, int Nr);
int ODEsolver();
int aorta(double* Param, int Nparam);

void runModel() {
	extern int Ncycles, Nmu, Nnu, Nr;

	initialgeometry(Nmu, Nnu, Nr);
	printf("Start ODEsolver with %i cycles ... ", Ncycles);
	double startTime = std::clock();
	ODEsolver();
	double endTime = std::clock();
	double programTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
	printf("run time = %f sec\n", programTime);
}

