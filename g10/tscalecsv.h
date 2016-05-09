#ifndef _TSCALECSV_H
#define _TSCALECSV_H

#include "derivs.h"
#include <stdio.h>

struct reaction {
	double (*k)(double, double[]);
	int n1, n2;
};

void tscaleoutput(FILE *fp, double *xp, double **yp, int kount);
double loadY(double *y, int idx);

#endif