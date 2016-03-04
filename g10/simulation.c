#include "simulation.h"

#define ZONES 1
#define NVAR 14*ZONES

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

void run(FILE *fp, double x2, bool print) {	
	dxsav = x2 / (KMAX * 1.2);
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_pX = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) vec_pX[i] = 0;
	
	double eps = 1e-5;
	double h1  = 1e-6;
	
	int nok = 0, nbad = 0;
	odeint(vec_pX, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
}