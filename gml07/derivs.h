#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>
#include <float.h>
#include "../constants.h"
// #include <cmath.h>
#include "../jacobian.h"

#define N_tot 300.0  // cm^-3
#define xC 1.4e-4
#define nC (N_tot*xC)
#define xO 3.2e-4
#define xSi 1.5e-5
#define nSi (N_tot*xSi)
#define xHe 0.1
#define xH_tot (1.0 - xC - xO - xSi - xHe)
#define N_H_tot (N_tot*xH_tot)

#define TEMP_INIT 50.0 // K
#define GRAIN_TEMP_INIT 50.0 // K

#define CDENSITY_H_tot 8e21 // cm^-2
#define SELF_SHIELDING (exp(-2e-21 * CDENSITY_H_tot) * 3.3e-11*1.7) // s^-1

#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1 try  try 1.8e-17

#define STEP_TIME 1 // s

double dnH2(double vec_pHx[]);
double dnH_p(double vec_pHx[]);

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getTemp();
void setTemp(int t);

double getGrainTemp();
void setGrainTemp(int t);

double getnH(double pH2, double pH_p);

double getxH2(double pH2);
double getxH_p(double pH_p);
double getxH(double pH2, double pH_p);

double getpH(double pH2, double pH_p);
double getpH2(double pH2);
double getpH_p(double pH_p);

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