#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>

#define nC 1.4e-4  // cm^-3
#define nO 3.2e-4  // cm^-3
#define nSi 1.5e-5 // cm^-3
#define nH_tot 91  // cm^-3

#define TEMP 600 // K
#define GRAIN_TEMP 50 // K

#define SELF_SHIELDING 1e-17*STEP_TIME // s^-1
#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1

#define STEP_TIME 1 // s

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getnH(double nH2, double nH_p);
double getne(double nH2, double nH_p);

double k1(double temp, double temp_grain);
double k2(double temp, double nH);
double k3(double temp, double nH);
double dk3ndH2(double temp, double nH2);
double k6(double temp);
double k7(double temp);
double k8(double temp);

#endif