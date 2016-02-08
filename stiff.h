#ifndef _STIFF_H_
#define _STIFF_H_

#include <math.h>
#include "nrutil.h"

void stiff(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], double *hdid, double *hnext, void (*derivs)(double, int, double[], double[]));
			
#endif