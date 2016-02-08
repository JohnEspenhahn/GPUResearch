#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>

#define TEMP 600 // K
#define GRAIN_TEMP 50 // K

#define STEP_TIME 1 // s

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getnH(double nH2, double nH_p);
double getne(double nH2, double nH_p);

double k1(double temp, double temp_grain);
double k2(double temp, double nH);
double k3(double temp, double nH);
double k6(double temp);
double k7(double temp);
double k8(double temp);

#endif