#include "odeint.h"
#include "main.h"
#include <stdio.h>

#define NVAR 2
#define END_STEPS 1000000000

#define nC 1.4e-4  // cm^-3
#define nO 3.2e-4  // cm^-3
#define nSi 1.5e-5 // cm^-3
#define nH_tot 91  // cm^-3

#define TEMP 1000 // K

double *vec_nHx;

double getnH2(int zone) { return vec_nHx[2*zone-1]; }
double getnH_p(int zone) { return vec_nHx[2*zone]; }
double getnH(int zone) { return nH_tot - getnH_p(zone) - getnH2(zone); }
double getne(int zone) { return getnH_p(zone) + nC + nSi; }

// t, nH2 are current values. ptr_dnH2dt is last slope. update ptr_dnH2dt
void dnH2dt(int zone, double t, double nH2, double *ptr_dnH2dt) {
	double nH = getnH(zone);
	
	*ptr_dnH2dt = k1(TEMP)*nH*nH - k2(TEMP,nH,nH2)*nH2*nH - k3(TEMP,nH,nH2)*nH2*nH2 - SELF_SHIELDING*nH2;
}

void dnH_pdt(int zone, double t, double nH_p, double *ptr_dnH_pdt) {
	double nH = getnH(zone), 
			ne = getne(zone);
	
	*ptr_dnH_pdt = COSMIC_RAY_RATE*nH + k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP)*nH_p*ne;
}

void dnHxdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i++) {
		if (IS_ODD(i)) { // if odd, nH2
			dnH2dt(i, t, vec_nHx[i], &vec_dnHxdt[i]);
		} else {
			dnH_pdt(i/2, t, vec_nHx[i], &vec_dnHxdt[i]);
		}
	}
}

int main() {
	vec_nHx = vector(1, NVAR);
	for (int i = 0; i < NVAR; i++) { vec_nHx[i] = 0; }
	
	double *vec_dnHxdt = vector(1, NVAR);
	dnHxdt(0, NVAR, vec_nHx, vec_dnHxdt);
	
	double x1 = 0, x2 = END_STEPS;
	double eps = 0.001;
	double h1 = 0.005;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, x1, x2, eps, h1, 0, &nok, &nbad, &dnHxdt);
	
	printf("integral(dnH2dt 0:%e) = %G\nintegral(dnH_pdt, 0:%e) = %G\n", (double) END_STEPS, vec_nHx[1], (double) END_STEPS, vec_nHx[2]);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n", nbad);
	
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