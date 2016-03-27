#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>
#include <float.h>
#include "../constants.h"
// #include <cmath.h>

#define N_TOT 300.0  // cm^-3
#define xC 1.4e-4
#define nC (N_TOT*xC)
#define xO 3.2e-4
#define xSi 1.5e-5
#define nSi (N_TOT*xSi)
#define xHe 0.1
#define xH_tot (1.0 - xC - xO - xSi - xHe)
#define N_H_tot (N_TOT*xH_tot)

#define TEMP_INIT 50 // K
#define GRAIN_TEMP_INIT 50 // K

#define SELF_SHIELDING 1e-17*STEP_TIME // s^-1
#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1 try  try 1.8e-17

#define STEP_TIME 1 // s

double dpH2(double vec_pHx[]);
double dpH_p(double vec_pHx[]);

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getTemp();
void setTemp(int t);

double getGrainTemp();
void setGrainTemp(int t);

double getnH(double pH2, double pH_p);
double getnH2(double pH2);
double getnH_p(double pH_p);

double getxH2(double pH2);
double getxH_p(double pH_p);
double getxH(double pH2, double pH_p);

double getpH(double pH2, double pH_p);

double getne(double pH_p);

double k1(double temp, double temp_grain);
double k2(double temp, double pH2, double pH_p);
double dk2dh2(double temp, double pH2, double pH_p);
double dk2dh_p(double temp, double pH2, double pH_p);
double k3(double temp, double pH2, double pH_p);
double dk3dh2(double temp, double pH2, double pH_p);
double dk3dh_p(double temp, double pH2, double pH_p);
double k6(double temp);
double k7(double temp);
double k8(double temp, double pH_p);
double dk8dph_p(double temp, double pH_p);

#endif