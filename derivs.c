#include "derivs.h"

double getnH(double nH2, double nH_p) { return nH_tot - nH_p - nH2; }
double getne(double nH2, double nH_p) { return nH_p + nC + nSi; }

// run a step of the derivatives
void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i += 2) {
		double nH2 = vec_nHx[i];
		double nH_p = vec_nHx[i+1];
		double nH = getnH(nH2, nH_p); 
		double ne = getne(nH2, nH_p);
		
		// nH2
		vec_dnHxdt[i] = k1(TEMP, GRAIN_TEMP)*nH*nH - k2(TEMP,nH)*nH2*nH - k3(TEMP,nH2)*nH2*nH2; // - SELF_SHIELDING*nH2;
		
		// nH_p
		vec_dnHxdt[i+1] = k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP)*nH_p*ne + COSMIC_RAY_RATE*nH;
	}
}

void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar)
{
	for (int i = 1; i <= nvar; i += 2) {
		dfdx[i] = 0.0;
		dfdx[i+1] = 0.0;
		
		double nH2 = vec_nHx[i];
		double nH_p = vec_nHx[i+1];
		double nH = getnH(nH2, nH_p);
		double ne = getne(nH2, nH_p);
		
		// nH2
		dfdy[i][1] = -k2(TEMP,nH)*nH - 2*k3(TEMP,nH)*nH2 - dk3ndH2(TEMP,nH2)*nH2*nH2; // need dk3ndH2 if k3 is actually function of nH2
		dfdy[i][2] = 0;
		
		// nH_p
		dfdy[i+1][1] = 0;
		dfdy[i+1][2] = -k7(TEMP)*ne - k8(TEMP)*ne;
	}
}

double k1(double temp, double temp_grain) { 
	double fa = 1/(1 + 1e4*exp(-600/temp_grain));
	double t2 = temp / 100;
	double t_grain2 = temp_grain / 100;
	
	return STEP_TIME * 3e-17 * ((sqrt(t2) * fa) / (1+0.4*sqrt(t2+t_grain2) + 0.2*t2 + 0.08*t2*t2));
}
double k2(double temp, double nH) {
	double kH = 1.2e-9*exp(-5.24e4/temp);
	double kL = 1.12e-10*exp(-7.035e4/temp);
	
	double log_temp = log(temp/1.0e4);
	double ncr = exp(4.0 - 0.416*log_temp - 0.327*log_temp*log_temp);
	
	return STEP_TIME * kH*pow(kH/kL, -ncr/(ncr+nH));
}
double k3(double temp, double nH2) { 
	double kH = 1.3e-9*exp(-5.33e4/temp);
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	double log_temp = log(temp/1.0e4);
	double ncr = exp(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	return STEP_TIME * kH*pow(kH/kL, -ncr/(ncr+nH2));
}
double dk3ndH2(double temp, double nH2) {
	double kH = 1.3e-9*exp(-5.33e4/temp);
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	double log_temp = log(temp/1.0e4);
	double ncr = exp(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	return STEP_TIME * k3(temp, nH2) * ((ncr*log(kH/kL)) / ((ncr+nH2)*(ncr+nH2)));
}
double k6(double temp) {
	return STEP_TIME * 5.466e-9*1.07*sqrt(temp/1.0e4)*exp(-13.6/(8.6173e-5*temp));
}
double k7(double temp) {
	return STEP_TIME * 2.54e-13*pow(temp/1.0e4, -0.8163 - 0.0208*log(temp/1.0e4));
}
double k8(double temp) {
	return STEP_TIME * 3.5e-12*pow(temp/300.0, -0.75);
}