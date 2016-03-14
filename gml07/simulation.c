#include "../simulation.h"
#include "derivs.h"

#define ZONES 1
#define NVAR 2*ZONES

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

void run(FILE *fp, double x2, bool print) {	
	dxsav = x2 / (KMAX * 1.2);
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_pHx = vector(1, NVAR);
	vec_pHx[1] = 5;
	vec_pHx[2] = 5;
	// for (int i = 1; i <= NVAR; i++) { vec_pHx[i] = 0; }
	
	double eps = 1e-5;
	double h1  = 1e-6;
	
	int nok = 0, nbad = 0;
	odeint(vec_pHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	double pH2 = vec_pHx[1];
	double pH_p = vec_pHx[2];
	double pH = getpH(pH2, pH_p);
	double ne = getne(pH_p);
	
	if (print) {
		printf("s     = %G\n\n", (double) x2);
		printf("%% pH2 = %G\n", getxH2(pH2) / (double) xH_tot);
		printf("OK calls %d\n", nok);
		printf("bad calls %d\n\n---------------------\n", nbad);
	}
	
	fprintf(fp, "x (Myr),%% H2,Temperature,Grain Temperature,%% H_p\n");	
	for (int i = 1; i < KMAX; i++) {
		// int i = KMAX - 1;
		double myr = xp[i]/(60*60*24*365*1e6);
		double perc_h2 = getxH2(yp[1][i]) / (double) xH_tot;
		double perc_hp = getxH_p(yp[2][i]) / (double) xH_tot;
		fprintf(fp, "%G,%G,%G,%G,%G\n", myr, perc_h2, getTemp(), getGrainTemp(), perc_hp);
	}
}