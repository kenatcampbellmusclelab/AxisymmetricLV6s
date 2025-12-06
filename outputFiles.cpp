// outputFiles.cpp
// TWS, November 2019, July 2025
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

void outputFiles(int run) {
	extern int Nt, Nnu, Nr, Nstore, Nparam, Nopt;
	extern double muin0, muout0, a0, nu_up, Ls0, Lsmax, Lsw, Tc;
	extern double rin0_la, rout0_la, rin0_rv, rout0_rv, rin0_ra, rout0_ra;
	extern double* tvec, * a1vec, * a2vec, * a3vec, * a4vec, * a5vec, * a6vec;
	extern double * Param;
	extern double** store, ** Y;
	extern std::string ParamName[112];	//note: zero-based arrays, include extra element


	int i, k, Nt1 = Nt - 1;
	char fname[80];
	FILE* ofp;
	fstream file1;

	// Write Matlab output file
	sprintf(fname, "outputFileA%03i.m", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "t = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", tvec[i]);
		else fprintf(ofp, "%g ]; \n", tvec[i]);
	}
	fprintf(ofp, "a1 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a1vec[i]);
		else fprintf(ofp, "%g ]; \n", a1vec[i]);
	}
	fprintf(ofp, "a2 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a2vec[i]);
		else fprintf(ofp, "%g ]; \n", a2vec[i]);
	}
	fprintf(ofp, "a3 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a3vec[i]);
		else fprintf(ofp, "%g ]; \n", a3vec[i]);
	}
	fprintf(ofp, "a4 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a4vec[i]);
		else fprintf(ofp, "%g ]; \n", a4vec[i]);
	}
	fprintf(ofp, "a5 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a5vec[i]);
		else fprintf(ofp, "%g ]; \n", a5vec[i]);
	}
	fprintf(ofp, "a6 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a6vec[i]);
		else fprintf(ofp, "%g ]; \n", a6vec[i]);
	}
	fclose(ofp);

	sprintf(fname, "outputFileStore%03i.m", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "store = [ ");
	for (k = 0; k < Nstore; k++) {
		for (i = 0; i < Nt1; i++) {
			if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
			if (i < Nt1 - 1) fprintf(ofp, "%g, ", store[k][i]);
			else if (k < Nstore - 1) fprintf(ofp, "%g; ...\n", store[k][i]);
			else fprintf(ofp, "%g ]; \n", store[k][i]);
		}
	}
	fclose(ofp);

	//Write combined text file - May 2022
	sprintf(fname, "outputFile_All%03i.txt", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "t a1 a2 a3 a4 a5 a6");
	fprintf(ofp, " Vtot Vlv Plv Vla Pla Vrv Prv Vra Pra null Psa Psp Psv Ppa Ppp Ppv q_la q_lv q_ra q_rv null q_sp q_sv null q_pp q_pv Zmiv Zaov");
	fprintf(ofp, " At     At_la  null  null null Ztcv Zpuv b-a_length dt(length) P_echo Q_trans P_ao_0 Q_lv_0");
	fprintf(ofp, " vel_lv vel_la vel_rv vel_ra LV_diam_int LV_diam_ext Vaorta\n");
	for (i = 0; i < Nt1; i++) {
		fprintf(ofp, "%g %g %g %g %g %g %g", tvec[i], a1vec[i], a2vec[i], a3vec[i], a4vec[i], a5vec[i], a6vec[i]);
		for (k = 0; k < Nstore; k++) fprintf(ofp, " %g", store[k][i]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}
