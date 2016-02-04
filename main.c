#include "odeint.h"
#include "derivs.h"
#include <stdio.h>

#define ZONES 1
#define NVAR 2*ZONES
#define END_STEPS 4000000
#define KMAX 255

FILE *fp;

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav = END_STEPS/KMAX;
int main() {
	fp = fopen("out.csv", "w");
	
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) { vec_nHx[i] = 0; }
	
	double x1 = 0, x2 = END_STEPS;
	double eps = 1e-5;
	double h1  = 1e-7;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, x1, x2, eps, h1, 0, &nok, &nbad, &dnHxdt, &bsstep);
	
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
	
	fprintf(fp, "x,dnH2dt,dnH_pdt\n");
	for (int i = 1; i < KMAX; i++) {
		fprintf(fp, "%G,%G,%G\n", xp[i], yp[1][i], yp[2][i]);
	}
	
	fclose(fp);
	
	return 0;
}