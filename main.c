#include "odeint.h"
#include "derivs.h"
#include "stiff.h"
#include <stdio.h>
#include <stdlib.h>

#define ZONES 1
#define NVAR 2*ZONES

#define KMAX 255

int kmax = KMAX, kount = 0;
double *xp, **yp, dxsav;
int main() {
	void run(int x2);
	
	char *p, s[100];
	int x2;
	
	do {
		printf("Number of seconds to run (0 to quit): ");
		p = 0;

		while (fgets(s, sizeof(s), stdin)) {
			x2 = strtol(s, &p, 10);
			if (p == s || *p != '\n') {
				printf("Please enter an integer: ");
			} else break;
		}
		
		if (x2 > 0) run(x2);
		else break;
	} while (1);
	
	return 0;
}

void run(int x2) {
	FILE *fp = fopen("out.csv", "w");
	
	dxsav = x2 / KMAX;
	xp = vector(1, KMAX);
	yp = matrix(1, NVAR, 1, KMAX);
	
	double *vec_nHx = vector(1, NVAR);
	for (int i = 1; i <= NVAR; i++) { vec_nHx[i] = 0; }
	
	double eps = 1e-4;
	double h1  = 1e-4;
	
	int nok = 0, nbad = 0;
	odeint(vec_nHx, NVAR, 0, x2, eps, h1, 0, &nok, &nbad, &derivs, &stiff);
	
	double nH2 = vec_nHx[1];
	double nH_p = vec_nHx[2];
	double nH = getnH(nH_p, nH2);
	double ne = getne(nH_p, nH2);
	
	printf("s    = %d\n\n", x2);
	printf("nH2  = %G\n", nH2);
	printf("nH_p = %G\n", nH_p);
	printf("nH   = %G\n", nH);
	printf("ne   = %G\n\n", ne);
	printf("OK calls %d\n", nok);
	printf("bad calls %d\n\n---------------------\n", nbad);
	
	fprintf(fp, "x,H2,H_p\n");
	for (int i = 1; i < KMAX; i++) {
		fprintf(fp, "%G,%G,%G\n", xp[i], yp[1][i], yp[2][i]);
	}
	
	fclose(fp);
}