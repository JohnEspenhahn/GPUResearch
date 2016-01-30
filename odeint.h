#ifndef _MMID_H_
#define _MMID_H_

#include <math.h>
#include "nrutil.h"

void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[], void (*derivs)(double, int, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
void bsstep(double y[], double dydx[], int nvar, double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, void (*derivs)(double, int, double[], double[]));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, void (*derivs)(double, int, double[], double[]));

#endif