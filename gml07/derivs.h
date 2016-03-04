#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>
#include <float.h>
#include "constants.h"
// #include <cmath.h>

#define N_TOT 300  // cm^-3
#define xC 1.4e-4
#define xO 3.2e-4
#define xSi 1.5e-5
#define xHe 0.1
#define xH_tot 1.0 - xC - xO - xSi - xHe

#define TEMP_INIT 60 // K
#define GRAIN_TEMP_INIT 45 // K

#define SELF_SHIELDING 1e-17*STEP_TIME // s^-1
#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1 try  try 1.8e-17

#define STEP_TIME 1 // s

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getnH_tot();
void setnH_tot(int t);

double getTemp();
void setTemp(int t);

double getGrainTemp();
void setGrainTemp(int t);

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