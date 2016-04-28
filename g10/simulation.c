#include "../simulation.h"
#include "derivs.h"

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;

void run(FILE *fp, double x2, bool print) {	
	dxsav = x2 / (KMAX * 1.2);
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) vec_nHx[i] = 1e-6;
	
	double eps = 1e-2;
	double h1  = 4;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	printf("%% nH2 = %G\n", getxH2(vec_nHx) / xH_tot);
	printf("%% nH  = %G\n", getxH(vec_nHx) / xH_tot);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n\n---------------------\n", nbad);
	
	fprintf(fp, "x (Myr),nH_p, nH2, nHe_p, nC_p, nO_p, nOH, nH2O, nCO , nC2, nO2 , nHCO_p, nCH, nCH2, nCH3_p, nH_m, nH2_p, nH3_p, nCH_p, nCH2_p, nOH_p, nH2O_p, nH3O_p, nCO_p, nHOC_p, nO_m, nC_m, nO2_p\n");
	for (int i = 1; i < KMAX; i++) {
		// int i = KMAX - 1;
		double myr = xp[i]/(60*60*24*365*1e6);
			 
		fprintf(fp, "%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G,%G\n"
		, myr, yp[1][i], yp[2][i], yp[3][i]
			 , yp[4][i], yp[5][i], yp[6][i]
			 , yp[7][i], yp[8][i], yp[9][i]
			 , yp[10][i], yp[11][i], yp[12][i]
			 , yp[13][i], yp[14][i], yp[15][i]
			 , yp[16][i], yp[17][i], yp[18][i]
			 , yp[19][i], yp[20][i], yp[21][i]
			 , yp[22][i], yp[23][i], yp[24][i]
			 , yp[25][i], yp[26][i], yp[27][i]);
	}
	
	free_vector(vec_nHx, 1, NVAR);
	free_matrix(yp, 1, NVAR, 1, KMAX);
	free_vector(xp, 1, KMAX);
}