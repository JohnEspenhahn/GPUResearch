#ifndef _MAIN_H_
#define _MAIN_H_

#include <math.h>

#define STEP_TIME 36000 // s

#define SELF_SHIELDING 1e-17*STEP_TIME // s^-1
#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1

void dnH2dt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void dnH_pdt(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);

double k1(double temp) { 
	return STEP_TIME * 6.6e-19 * sqrt(temp); 
}
double k2(double temp, double nH, double nH2) {
	double kH = 1.2e-9*exp(-5.24e4/temp);
	double kL = 1.12e-10*exp(-7.035e4/temp);
	
	double log_temp = log(temp/1.0e4);
	double ncr = pow(4.0 - 0.416*log_temp - 0.327*log_temp*log_temp, 10);
	
	return STEP_TIME * kH*pow(kH/kL, -1.0/(1+nH/ncr));
}
double k3(double temp, double nH, double nH2) { 
	double kH = 1.3e-9*exp(-5.33e4/temp);
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	double log_temp = log(temp/1.0e4);
	double ncr = pow(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp, 10);
	
	return STEP_TIME * kH*pow(kH/kL, -1.0/(1+nH/ncr));
}
double k6(double temp) {
	double logT = log(temp);
	
	return STEP_TIME * exp(-32.71396786 
				+ 13.536556*logT 
				- 5.73932875*logT*logT 
				+ 1.56315498*logT*logT*logT 
				- 0.2877056*logT*logT*logT*logT
				+ 3.48255977e-2*logT*logT*logT*logT*logT
				- 2.63197617e-3*logT*logT*logT*logT*logT*logT
				+ 1.11954395e-4*logT*logT*logT*logT*logT*logT*logT
				- 2.03914985e-6*logT*logT*logT*logT*logT*logT*logT*logT);
}
double k7(double temp) {
	return STEP_TIME * 2.54e-13*pow(temp/1.0e4, -0.8163 - 0.0208*log(temp/1.0e4));
}
double k8(double temp) {
	return STEP_TIME * 3.5e-12*pow(temp/300.0, -0.75);
}

#endif