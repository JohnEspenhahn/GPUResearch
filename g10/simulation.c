#include "../simulation.h"
#include "derivs.h"
#include "tscalecsv.h"

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

FILE *csv_fp;

void run(FILE *fp, double x2, double eps, double h1, bool print) {	
	dxsav = x2 / KMAX;
	xp = vector(1, KMAX);
	yp = matrix(1, KMAX, 1, NVAR);
	csv_fp = fp;
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) vec_nHx[i] = 1e-6;
	
	// double eps = 1e8;
	// double h1  = 1e3;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	printf("%% nH2 = %G\n", getxH2(vec_nHx) / xH_tot);
	printf("%% nH  = %G\n", getxH(vec_nHx) / xH_tot);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n\n---------------------\n", nbad);
	
	tscaleoutput(fp, xp, yp, kount);
	
	free_vector(vec_nHx, 1, NVAR);
	free_matrix(yp, 1, KMAX, 1, NVAR);
	free_vector(xp, 1, KMAX);
}