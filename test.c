#include "derivs.h"
#include "nrutil.h"
#include <stdio.h>

#define TEMP 700 // K
#define GRAIN_TEMP 75 // K

#define KMAX 0

FILE *fp;

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav = 100;

void outputRates() {
	fp = fopen("out.csv", "w");
	fprintf(fp, "dnH2dt,nH,nH2,dnH_pdt,ne,nH_p\n");
	
	double *vec_nHx = vector(1,2), *vec_dnHxdt = vector(1,2);
	vec_nHx[1] = vec_nHx[2] = vec_dnHxdt[1] = vec_dnHxdt[2] = 0;
	
	double targ = 100;
	double steps = 255;
	double step_amnt = targ / steps;
	
	for (int i = 0; i < steps; i++) {
		dnHxdt(0, 2, vec_nHx, vec_dnHxdt);
		
		vec_nHx[1] += step_amnt;
		vec_nHx[2] += step_amnt;
	}
	
	fclose(fp);
}

int main() {
	double nH[11] = { 
					-16.451,
					-17.2353,
					-18.0196,
					-18.8039,
					-19.5882,
					-20.3725,
					-21.1569,
					-21.9412,
					-22.7255,
					-23.5098,
					-24.2941
				};
	double nH2[11] = {
					53.7255,
					54.1176,
					54.5098,
					54.902,
					55.2941,
					55.6863,
					56.0784,
					56.4706,
					56.8627,
					57.2549,
					57.6471
				};
				
	for (int i = 0; i < 11; i++) {
		printf("k1 %G\n", k1(TEMP, GRAIN_TEMP)*nH[i]*nH[i]);
		printf("k2 %G\n", -k2(TEMP,nH[i],nH2[i])*nH2[i]*nH[i]);
		printf("k3 %G\n\n", -k3(TEMP,nH[i],nH2[i])*nH2[i]*nH2[i]);
	}
	
	return 0;
}