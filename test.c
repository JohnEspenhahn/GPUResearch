#include "test.h"
#include "gml07/derivs.h"
#include "odeint.h"
#include "stiff.h"
#include <stdio.h>
#include "simulation.h"

/*
void plotK1(FILE *fp, int steps) {
	double start_temp = 10, end_temp = 200;
	double start_grain_temp = 10, end_grain_temp = 60;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_grain_temp; j < end_grain_temp; j += (end_grain_temp - start_grain_temp) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, k1(i,j));
		}
	}
}

void plotK2(FILE *fp, int steps) {
	double start_temp = 10, end_temp = 200;
	double start_nH = 0, end_nH = 200;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_nH; j < end_nH; j += (end_nH - start_nH) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, k2(i,j));
		}
	}
}

void plotK3(FILE *fp, int steps) {
	double start_temp = 10, end_temp = 200;
	double start_nH2 = 0, end_nH2 = 200;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_nH2; j < end_nH2; j += (end_nH2 - start_nH2) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, k3(i,j));
		}
	}
}

void plotdk3ndH2(FILE *fp, int steps) {
	double start_temp = 10, end_temp = 200;
	double start_nH2 = 0, end_nH2 = 200;
	for (double i = start_temp; i < end_temp; i += (end_temp - start_temp) / steps) {
		for (double j = start_nH2; j < end_nH2; j += (end_nH2 - start_nH2) / steps) {
			fprintf(fp, "%G,%G,%G\n", i, j, dk3ndH2(i,j));
		}
	}
}

void plotK6(FILE *fp, int steps) { 
	for (double i = 200; i < 1000; i += (1000.0-200.0)/steps) {
		fprintf(fp, "%G,%G\n", i, k6(i));
	}
}

void plotH_p(FILE *fp, int steps) {
	double nH = 400;
	for (double t = 50; t < 600; t += 10) {
		double nH_p = 1e-4;
		double ne = 1e-4;
		fprintf(fp, "%G,%G,%G\n", nH, t, k6(t)*nH*ne - k7(t)*nH_p*ne - k8(t)*nH_p*ne);
	}
}

void plotH2(FILE *fp, int steps) {
	double nH = 400;
	for (double t = 50; t < 600; t += 5) {
		for (double gt = 20; gt < 100; gt += 5) {
			double nH2 = 100;
			fprintf(fp, "%G,%G,%G,%G\n", nH, t, gt, k1(t, gt)*nH*nH - k2(t,nH)*nH2*nH - k3(t,nH2)*nH2*nH2);
		}
	}
}
*/

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
			
			fprintf(fp, "%G,%G,%G,%G\n", dfdy[1][1], dfdy[1][2], dfdy[2][1], dfdy[2][2]);
		}
	}
}

void plotRuns(FILE *fp, int steps) {
	int init = 75, max = 200;
	for (int j = init; j < max; j += (max - init) / steps) {
		setTemp(j);
		run(fp, 5e14, 0);
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	// fprintf(fp, "temp,grain_temp,k1,k2,k3,k6,k7,k8\n");
	
	int steps = 75;
	plotJacobn(fp, steps);
	
	// dk3ndH2(73.3333, 200);
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}