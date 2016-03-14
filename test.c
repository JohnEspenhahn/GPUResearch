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
			fprintf(fp, "%G,%G,%G,%G\n", nH, t, gt, k1(TEMP, GRAIN_TEMP)*pH*nH - k2(TEMP,pH2,pH_p)*pH2*nH - k3(TEMP,pH2,pH_p)*pH2*nH2 + SELF_SHIELDING*pH2);
		}
	}
}
*/

void plotDerivs(FILE *fp) {
	double *nHx = vector(1,2);
	double *dnHx = vector(1,2);
	
	for (double pH2 = 0; pH2 < 90; pH2 += 2) {
		for (double pH_p = 0; pH_p < 90; pH_p += 2) {
			nHx[1] = pH2;
			nHx[2] = pH_p;
			
			derivs(0, 2, nHx, dnHx);
			fprintf(fp, "%G,%G,%G,%G\n", pH2, pH_p, dnHx[1], dnHx[2]);
		}
	}
}

void plotJacobn(FILE *fp) {
	double *nHx = vector(1,2);	
	double *dfdx = vector(1,2);
	double **dfdy = matrix(1,2,1,2);
	
	for (double pH2 = 0; pH2 < 90; pH2 += 2) {
		for (double pH_p = 0; pH_p < 90; pH_p += 2) {
			nHx[1] = pH2;
			nHx[2] = pH_p;
			
			jacobn(0, nHx, dfdx, dfdy, 2);
			fprintf(fp, "%G,%G,%G,%G,%G,%G\n", pH2, pH_p, dfdy[1][1], dfdy[1][2], dfdy[2][1], dfdy[2][2]);
		}
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	
	int steps = 75;
	plotDerivs(fp);
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}