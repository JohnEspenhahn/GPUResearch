#include "test.h"
#include "derivs.h"
#include "odeint.h"
#include "stiff.h"
#include <stdio.h>
#include "simulation.h"

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
	double max_H = 200;
	double max_t = 700;
	for (double nH = 10; nH < max_H; nH += max_H / 75.0) { // percent hydrogen
		for (double t = 50; t < max_t; t += max_t / 75.0) {
			double nH_p = max_H - nH;
			double nH2 = max_H / 2.0;
			fprintf(fp, "%G,%G,%G\n", nH, t, k1(t, getGrainTemp())*nH*nH - k2(t,nH)*nH2*nH - k3(t,nH2)*nH2*nH2);
		}
	}
}

void plotJacobn(FILE *fp, int steps) {
	double *vec_nHx = vector(1,2);
	vec_nHx[1] = 1e1;
	vec_nHx[2] = 1e-3;
	
	double *dfdx = vector(1,2);
	double **dfdy = matrix(1,2,1,2);
	
	for (int temp = 1; temp < 700; temp += 20) {
		for (int grain = 1; grain < 100; grain += 10) {
			setTemp(temp);
			setGrainTemp(grain);
			jacobn(0, vec_nHx, dfdx, dfdy, 2);
			
			fprintf(fp, "%G,%G\n", dfdy[1][1], dfdy[2][2]);
		}
	}
}

void plotRuns(FILE *fp, int _) {
	int init = 70, max = 1000;
	for (int j = init; j < max; j += (max - init) / 30) {
		for (int i = 10; i < 100; i += 100 / 25) {
			setTemp(j);
			setGrainTemp(i);
			run(fp, 5e14);
		}
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	// fprintf(fp, "temp,grain_temp,k1,k2,k3,k6,k7,k8\n");
	
	int steps = 75;
	plotRuns(fp, steps);
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}