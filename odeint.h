#ifndef _MMID_H_
#define _MMID_H_

#include <math.h>
#include "nrutil.h"

void mmid(float y[], float dydx[], int nvar, float xs, float htot, int nstep, float yout[], void (*derivs)(float, int, float[], float[]));
void pzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);
void bsstep(float y[], float dydx[], int nvar, float *xx, float htry, float eps, float yscal[], float *hdid, float *hnext, void (*derivs)(float, int, float[], float[]));
void odeint(float ystart[], int nvar, float x1, float x2, float eps, float h1, float hmin, int *nok, int *nbad, void (*derivs)(float, int, float[], float[]));

#endif