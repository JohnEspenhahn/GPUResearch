#include "odeint.h"
#include "derivs.h"
#include <stdio.h>

#define ZONES 1
#define NVAR 2*ZONES
#define END_STEPS 4000000
#define KMAX 0

FILE *fp;

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav = 100;
int main() {
	fp = fopen("out.csv", "w");
	fprintf(fp, "dnH2dt,nH,nH2,dnH_pdt,ne,nH_p\n");
	
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 0; i < NVAR; i++) { vec_nHx[i] = 0; }
	
	double x1 = 0, x2 = END_STEPS;
	double eps = 1e-4;
	double h1  = 1e-6;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, x1, x2, eps, h1, 0, &nok, &nbad, &dnHxdt);
	
	double nH2 = vec_nHx[1];
	double nH_p = vec_nHx[2];
	double nH = getnH(nH_p, nH2);
	double ne = getne(nH_p, nH2);
	
	printf("nH2  = %G\n", nH2);
	printf("nH_p = %G\n", nH_p);
	printf("nH   = %G\n", nH);
	printf("ne   = %G\n", ne);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n", nbad);
	
	fclose(fp);
	
	return 0;
}