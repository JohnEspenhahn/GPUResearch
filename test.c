#include "test.h"
#include "derivs.h"
#include "odeint.h"
#include "stiff.h"
#include <stdio.h>

#define KMAX 0

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav = 100;

void plotK1(FILE *fp, int steps) {
	double start_temp = 200, end_temp = 1000;
	double start_grain_temp = 10, end_grain_temp = 100;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_grain_temp; j < end_grain_temp; j += (end_grain_temp - start_grain_temp) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, k1(i,j));
		}
	}
}

void plotK2(FILE *fp, int steps) {
	double start_temp = 595, end_temp = 605;
	double start_nH = 20, end_nH = 100;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_nH; j < end_nH; j += (end_nH - start_nH) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, k2(i,j));
		}
	}
}

void plotK6(FILE *fp, int steps) { 
	for (double i = 200; i < 1000; i += (1000.0-200.0)/steps) {
		fprintf(fp, "%G,%G\n", i, k6(i));
	}
}

void plotH_p(FILE *fp, int steps) {
	double max_H = 100;
	for (double nH = 0; nH < max_H; nH += max_H / 75.0) { // percent hydrogen
		for (double t = 200; t < 1000; t += 900.0 / 75.0) {
			double nH_p = max_H - nH;
			double ne = getne(0, nH_p);
			fprintf(fp, "%G,%G,%G\n", nH, t, k6(t)*nH*ne - k7(t)*nH_p*ne - k8(t)*nH_p*ne);
		}
	}
}

void plotH2(FILE *fp, int steps) {
	double s_H = 0, e_H = 90;
	double s_H2 = 0, e_H2 = 90;
	for (double nH = s_H; nH < e_H; nH += (e_H - s_H) / steps) {
		for (double nH2 = s_H2; nH2 < e_H2; nH2 += (e_H2 - s_H2) / steps) {
			fprintf(fp, "%G,%G,%G\n", nH, nH2, k1(TEMP, GRAIN_TEMP)*nH*nH - k2(TEMP,nH)*nH2*nH - k3(TEMP,nH2)*nH2*nH2);
		}
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	// fprintf(fp, "temp,grain_temp,k1,k2,k3,k6,k7,k8\n");
	
	int steps = 75;
	plotH_p(fp, steps);
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}