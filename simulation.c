#include "simulation.h"

#define ZONES 1
#define NVAR 2*ZONES

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

void run(FILE *fp, double x2, bool print) {	
	dxsav = x2 / (KMAX * 1.2);
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	vec_nHx[1] = 0;
	vec_nHx[2] = 0;
	// for (int i = 1; i <= NVAR; i++) { vec_nHx[i] = 0; }
	
	double eps = 1e-5;
	double h1  = 1e-6;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	double nH2 = vec_nHx[1];
	double nH_p = vec_nHx[2];
	double nH = getnH(nH2, nH_p);
	double ne = getne(nH2, nH_p);
	
	if (print) {
		printf("s     = %G\n\n", (double) x2);
		printf("%% nH2 = %G\n", nH2 / (double) getnH_tot());
		printf("OK calls %d\n", nok);
		printf("bad calls %d\n\n---------------------\n", nbad);
	}
	
	fprintf(fp, "x (Myr),%% H2,H Density,Temperature,Grain Temperature,%% H_p\n");	
	for (int i = 1; i < KMAX; i++) {
		// int i = KMAX - 1;
		double myr = xp[i]/(60*60*24*365*1e6);
		double perc_h2 = yp[1][i] / (double) getnH_tot();
		double perc_hp = yp[2][i] / (double) getnH_tot();
		fprintf(fp, "%G,%G,%G,%G,%G,%G\n", myr, perc_h2, getnH_tot(), getTemp(), getGrainTemp(), perc_hp);
	}
}