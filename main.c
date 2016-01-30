#include "odeint.h"
#include "main.h"
#include <stdio.h>

#define NVAR 1
#define END_STEPS 100000

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

// t, nH2 are current values. dnH2dt is last slope. update dnH2dt
void dnH2dt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int zone = 1; zone <= nvar; zone++) {
		double nH = getnH(zone), 
		      nH2 = getnH2(zone);
		
		vec_dnHxdt[2*zone-1] = k1(TEMP)*nH*nH - k2(TEMP,nH,nH2)*nH2*nH - k3(TEMP,nH,nH2)*nH2*nH2 - SELF_SHIELDING*nH2;
	}
}

void dnH_pdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int zone = 1; zone <= nvar; zone++) {
		double nH = getnH(zone), 
		      nH_p = getnH_p(zone),
			  ne = getne(zone);
		
		vec_dnHxdt[2*zone] = COSMIC_RAY_RATE*nH + k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP)*nH_p*ne;
	}
}

void dnHxdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	dnH2dt(t, nvar, vec_nHx, vec_dnHxdt);
	// dnH_pdt(t, nvar, vec_nHx, vec_dnHxdt);
}

int main() {
	vec_nHx = vector(1, NVAR);
	for (int i = 0; i < NVAR; i++) { vec_nHx[i] = 0; }
	
	double *vec_dnHxdt = vector(1, NVAR);
	dnHxdt(0, NVAR, vec_nHx, vec_dnHxdt);
	
	double x1 = 0, x2 = END_STEPS;
	double eps = 0.01;
	double h1 = 0.005;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, x1, x2, eps, h1, 0, &nok, &nbad, &dnHxdt);
	
	printf("integral(dnH2dt 0:%e) = %G\nintegral(dnH_pdt, 0:%e) = %G", (double) END_STEPS, vec_nHx[1], (double) END_STEPS, vec_nHx[2]);
	
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