#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>

#define STEP_TIME 1 // s

void dnHxdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);

void dnH2dt(double t, double nH_p, double nH2, double *ptr_dnH2dt);
void dnH_pdt(double t, double nH_p, double nH2, double *ptr_dnH_pdt);

double getnH(double nH_p, double nH2);
double getne(double nH_p, double nH2);

double k1(double temp, double temp_grain);
double k2(double temp, double nH, double nH2);
double k3(double temp, double nH, double nH2);
double k6(double temp);
double k7(double temp);
double k8(double temp);

#endif