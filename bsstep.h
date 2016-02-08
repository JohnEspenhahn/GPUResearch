#ifndef _BSSTEP_H_
#define _BSSTEP_H_

#include <math.h>
#include "nrutil.h"

void bsstep(double y[], double dydx[], int nvar, double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, void (*derivs)(double, int, double[], double[]));

#endif