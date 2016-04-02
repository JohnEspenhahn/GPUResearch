#include "test.h"
#include "gml07/derivs.h"
#include "odeint.h"
#include "stiff.h"
#include <stdio.h>
#include "simulation.h"
#include "jacobian.h"

double end_pH2 = 1e-22;
double end_pH_p = 1e-25;
double end_temp = 199;

void plotK1(FILE *fp, int steps) {
	for (double temp = 1; temp <= 1001; temp += 100) {
		for (double grain_temp = 1; grain_temp <= 1001; grain_temp += 100) {
			fprintf(fp, "k1,%G,%G,,%G\n", temp, grain_temp, k1(temp,grain_temp));
		}
	}
}

void plotK2(FILE *fp, int steps) {
	for (double pH2 = 0; pH2 < end_pH2; pH2 += end_pH2 / steps) {
		for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p / steps) {
			fprintf(fp, "k2,%G,%G,%G,%G\n", TEMP_INIT, pH2, pH_p, k2(TEMP_INIT,pH2,pH_p));
		}
	}
}

void plotK3(FILE *fp, int steps) {
	for (double pH2 = 0; pH2 < end_pH2; pH2 += end_pH2 / steps) {
		for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p / steps) {
			fprintf(fp, "k3,%G,%G,%G,%G\n", TEMP_INIT, pH2, pH_p, k3(TEMP_INIT,pH2,pH_p));
		}
	}
}

void plotK6(FILE *fp, int steps) { 
	fprintf(fp, "k6,%G,,,%G\n", TEMP_INIT, k6(TEMP_INIT));
}

void plotK7(FILE *fp, int steps) { 
	fprintf(fp, "k7,%G,,,%G\n", TEMP_INIT, k7(TEMP_INIT));
}

void plotK8(FILE *fp, int steps) { 
	for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p / steps) {
		fprintf(fp, "k8,%G,,%G,%G\n", TEMP_INIT, pH_p, k8(TEMP_INIT,pH_p));
	}
}

void plotH_p(FILE *fp, int steps) {
	for (double temp = 1; temp < end_temp; temp += end_temp / steps) {
		for (double pH2 = 0; pH2 < end_pH2; pH2 += end_pH2 / steps) {
			for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p / steps) {
				double ne = getne(pH_p)
					 , pH = getpH(pH2,pH_p);
				fprintf(fp, "dH_p,%G,%G,%G,%G\n", temp, pH2, pH_p, k6(temp)*pH*ne - k7(temp)*pH_p*ne - k8(temp,pH_p)*pH_p*ne + COSMIC_RAY_RATE*pH);
			}
		}
	}
}

void plotH2(FILE *fp, int steps) {
	for (double temp = 1; temp < end_temp; temp += end_temp / steps) {
		for (double pH2 = 0; pH2 < end_pH2; pH2 += end_pH2 / steps) {
			for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p / steps) {
				double nH = getnH(pH2,pH_p);
				fprintf(fp, "dH2,%G,%G,%G,%G\n", temp, pH2, pH_p, k1(temp,temp)*getpH(pH2,pH_p)*nH - k2(temp,pH2,pH_p)*pH2*nH - k3(temp,pH2,pH_p)*pH2*getnH2(pH2) - SELF_SHIELDING*pH2);
			}
		}
	}
}

void plotDerivs(FILE *fp) {
	double *nHx = vector(1,2);
	double *dnHx = vector(1,2);
	
	for (double pH2 = 2; pH2 < 90; pH2 += 2) {
		for (double pH_p = 2; pH_p < 90; pH_p += 2) {
			nHx[1] = pH2;
			nHx[2] = pH_p;
			
			derivs(0, 2, nHx, dnHx);
			fprintf(fp, "%G,%G,%G,%G\n", pH2, pH_p, dnHx[1], dnHx[2]);
		}
	}
}

void compareJacobn(FILE *fp, int steps) {
	double *nHx = vector(1,2);	
	double *dfdx = vector(1,2);
	double **dfdy = matrix(1,2,1,2);
	double (*derivs[])(double[]) = { &dpH2, &dpH_p };
	
	for (double pH2 = 0; pH2 < end_pH2; pH2 += end_pH2/steps) {
		for (double pH_p = 0; pH_p < end_pH_p; pH_p += end_pH_p/steps) {
			nHx[1] = pH2;
			nHx[2] = pH_p;
			
			// fprintf(fp, "nH2, nH_p, nH = %G,%G,%G\n", getnH2(pH2), getnH_p(pH_p), getnH(pH2,pH_p));
			printf("pH2, pH_p, pH = %G,%G,%G\n", pH2, pH_p, getpH(pH2,pH_p));
			
			jacobn(0, nHx, dfdx, dfdy, 2);
			printf("%G,%G,%G,%G\n", dfdy[1][1], dfdy[1][2], dfdy[2][1], dfdy[2][2]);
			
			jacobian(derivs, 2, nHx, 1e-30, dfdy, 2);
			printf("%G,%G,%G,%G\n", dfdy[1][1], dfdy[1][2], dfdy[2][1], dfdy[2][2]);
		}
	}
}

void printArray(double* arr, int i, int j) {
	for (;i <= j;i++) {
		printf("%f\t", arr[i]);
	}
	printf("\n");
}

void testCopy() {
	double *v1 = vector(1,5), *v2 = vector(1,10);
	for (int i = 1; i <= 10; i++) v2[i] = i;
	
	copy_vector(v2+5, v1, 1, 5);
	
	printArray(v1, 1, 5);
	printArray(v2, 1, 10);
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	
	int steps = 5;	
	
	// printf("f sheild %G\n", fshield);
	printf("self sheild %G\n", SELF_SHIELDING);
	
	// plotK1(fp, steps);
	
	// plotH2(fp, steps);
	// plotH_p(fp, steps);
	
	/* plotK1(fp, steps);
	plotK2(fp, steps);
	plotK3(fp, steps);
	plotK6(fp, steps);
	plotK7(fp, steps);
	plotK8(fp, steps);	 */
	
	// compareJacobn(fp, steps);
	// testCopy();
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}