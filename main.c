#include "odeint.h"
#include "main.h"
#include <stdio.h>

#define ZONES 1
#define NVAR 2*ZONES
#define END_STEPS 1000
#define KMAX 0

#define nC 1.4e-4  // cm^-3
#define nO 3.2e-4  // cm^-3
#define nSi 1.5e-5 // cm^-3
#define nH_tot 91  // cm^-3

#define TEMP 700 // K
#define GRAIN_TEMP 75 // K

long H2_calls = 0, H_p_calls = 0;

FILE *fp;

double getnH(double nH_p, double nH2) { return nH_tot - nH_p - nH2; }
double getne(double nH_p) { return nH_p + nC + nSi; }

// t, nH2 are current values. ptr_dnH2dt is last slope. update ptr_dnH2dt
void dnH2dt(double t, double nH_p, double nH2, double *ptr_dnH2dt) {
	double nH = getnH(nH_p, nH2);
	
	H2_calls += 1;
	*ptr_dnH2dt = k1(TEMP, GRAIN_TEMP)*nH*nH - k2(TEMP,nH,nH2)*nH2*nH - k3(TEMP,nH,nH2)*nH2*nH2; // - SELF_SHIELDING*nH2;
	
	fprintf(fp, "%G,%G,%G,", (*ptr_dnH2dt), nH, nH2);
}

void dnH_pdt(double t, double nH_p, double nH2, double *ptr_dnH_pdt) {
	double nH = getnH(nH_p, nH2), 
			ne = getne(nH_p);
	
	H_p_calls += 1;
	*ptr_dnH_pdt = k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP)*nH_p*ne; // + COSMIC_RAY_RATE*nH;
	
	fprintf(fp, "%G,%G,%G\n", (*ptr_dnH_pdt), ne, nH_p);
}

void dnHxdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i++) {
		if (IS_ODD(i)) { // if odd, nH2
			double nH2 = vec_nHx[i];
			double nH_p = vec_nHx[i+1];
			dnH2dt(t, nH_p, nH2, &vec_dnHxdt[i]);
		} else {
			double nH_p = vec_nHx[i];
			double nH2 = vec_nHx[i-1];
			dnH_pdt(t, nH_p, nH2, &vec_dnHxdt[i]);
		}
	}
}

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav = 100;
int main() {
	fp = fopen("out.csv", "w");
	fprintf(fp, "dnH2,nH,nH2,dnH_p,ne,nH_p\n");
	
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 0; i < NVAR; i++) { vec_nHx[i] = 0; }
	
	double x1 = 0, x2 = END_STEPS;
	double eps = 0.001;
	double h1 = 0.005;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, x1, x2, eps, h1, 0, &nok, &nbad, &dnHxdt);
	
	printf("\nH2 calls %ld\n", H2_calls);
	printf("H_p calls %ld\n", H_p_calls);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n", nbad);
	
	fclose(fp);
	
	return 0;
}

/*
int main_xcube() {
	double *y, *dydx;
	y = vector(1,1);
	y[1] = 0;
	
	dydx = vector(1,1);
	xcube(0, y, dydx);
	
	int xs = 0, htot = 2, nstep = 10;
	mmid(y, dydx, 1, xs, htot, nstep, y, *xcube);
	
	printf("integral(3x^2, %i:%f:%i) = %f\n", xs, (double) htot/nstep, xs+htot, y[1]);
	
	return 0;
}
*/