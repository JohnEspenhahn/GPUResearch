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

void plotH_p(FILE *fp, int steps) {
	double s_H = 55, e_H = 65;
	double s_Hp = 1, e_Hp = 91;
	for (double nH = s_H; nH < e_H; nH += (e_H - s_H) / steps) {
		for (double nH_p = s_Hp; nH_p < e_Hp; nH_p += (e_Hp - s_Hp) / steps) {
			double ne = getne(0, nH_p);
			fprintf(fp, "%G,%G,%G\n", nH, nH_p, k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP)*nH_p*ne);
		}
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	// fprintf(fp, "temp,grain_temp,k1,k2,k3,k6,k7,k8\n");
	
	int steps = 50;
	plotH_p(fp, steps);
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}