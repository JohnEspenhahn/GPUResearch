#include "../simulation.h"
#include "derivs.h"

#define ZONES 1
#define NVAR 14*ZONES

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

void run(FILE *fp, double x2, bool print) {	
	dxsav = x2 / (KMAX * 1.2);
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) vec_nHx[i] = 1e-5;
	
	double eps = 1e-4;
	double h1  = 1e-4;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	printf("%% nH2 = %G\n", getxH2(vec_nHx) / xH_tot);
	printf("%% nH  = %G\n", getxH(vec_nHx) / xH_tot);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n\n---------------------\n", nbad);
	
	free_vector(vec_nHx, 1, NVAR);
	free_matrix(yp, 1, NVAR, 1, KMAX);
	free_vector(xp, 1, KMAX);
}