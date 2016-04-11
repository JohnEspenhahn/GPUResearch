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
	
	double *vec_nHx = vector(1, NVAR);
	vec_nHx[1] = 1e-3;
	vec_nHx[2] = 1e-3;
	// for (int i = 1; i <= NVAR; i++) { vec_nHx[i] = 0; }
	
	double eps = 1e-9;
	double h1  = 1e-5;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	double nH2 = vec_nHx[1];
	double nH_p = vec_nHx[2];
	
	if (print) {
		printf("s     = %G\n\n", (double) x2);
		printf("pH2, pH_p, pH = %G,%G,%G\n", getpH2(nH2), getpH_p(nH_p), getpH(nH2,nH_p));
		printf("nH2, nH_p, nH = %G,%G,%G\n", nH2, nH_p, getpH(nH2,nH_p));
		printf("sum pHs = %G\n", getpH2(nH2)+getpH_p(nH_p)+getpH(nH2,nH_p));
		printf("sum nHs = %G\n", nH2*mu_h2*M_h + nH_p*mu_h*M_h + getnH(nH2,nH_p)*mu_h*M_h);
		
		printf("%% pH2 = %G\n", getxH2(nH2) / xH_tot);
		printf("%% pHp = %G\n", getxH_p(nH_p) / xH_tot);
		printf("%% pH  = %G\n", getxH(nH2, nH_p) / xH_tot);
		printf("sum %% pHs = %G\n", (getxH2(nH2)+getxH_p(nH_p)+getxH(nH2, nH_p)) / xH_tot);
		printf("OK calls %d\n", nok);
		printf("bad calls %d\n\n---------------------\n", nbad);
	}
	
	fprintf(fp, "x (Myr),%% H2,Temperature,Grain Temperature,%% H_p\n");	
	for (int i = 1; i < KMAX; i++) {
		// int i = KMAX - 1;
		double myr = xp[i]/(60*60*24*365*1e6);
		double perc_h2 = getxH2(yp[1][i]);
		double perc_hp = getxH_p(yp[2][i]);
		fprintf(fp, "%G,%G,%G,%G,%G\n", myr, perc_h2, getTemp(), getGrainTemp(), perc_hp);
	}
	
	free_vector(vec_nHx, 1, NVAR);
	free_matrix(yp, 1, NVAR, 1, KMAX);
	free_vector(xp, 1, KMAX);
}