#include "derivs.h"

/* private */ int TEMP = TEMP_INIT, GRAIN_TEMP = GRAIN_TEMP_INIT;

double getTemp() { return TEMP; }
void setTemp(int t) { if (t > 0) TEMP = t; }

double getGrainTemp() { return GRAIN_TEMP; }
void setGrainTemp(int t) { if (t > 0) GRAIN_TEMP = t; }

double getnH(double nH2, double nH_p) { return N_H_tot - 2*nH2 - nH_p; }

double getpH2(double nH2) { return nH2 * (mu_h2*M_h); }
double getpH_p(double nH_p) { return nH_p * (mu_h*M_h); }
double getpH(double nH2, double nH_p) { return getnH(nH2, nH_p)*(mu_h*M_h); };

double getxH2(double nH2) { return (2*nH2) / N_tot; }
double getxH_p(double nH_p) { return nH_p / N_tot; }
double getxH(double nH2, double nH_p) { return xH_tot - getxH2(nH2) - getxH_p(nH_p); }

/** Derives the abundance of e- based on the mass densities of H2 and H+ */
double getxe(double nH_p) { 
	return getxH_p(nH_p) + xC + xSi; 
}

/** Derives the number density of e- based on the mass density of H+ */
double getne(double nH_p) {
	return getxe(nH_p) * N_tot;
}

double dnH2(double vec_nHx[]) {
	double nH2 = vec_nHx[1]
		 , nH_p = vec_nHx[2]
		 , nH = getnH(nH2, nH_p);
	
	return (k1(TEMP, GRAIN_TEMP)*nH*nH - k2(TEMP,nH2,nH_p)*nH2*nH - k3(TEMP,nH2,nH_p)*nH2*nH2 - SELF_SHIELDING*nH2);
}

double dnH_p(double vec_nHx[]) {
	double nH2  = vec_nHx[1]
	     , nH_p = vec_nHx[2]
		 , nH   = getnH(nH2, nH_p)
		 , ne   = getne(nH_p);
	
	return (k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP,ne)*nH_p*ne + COSMIC_RAY_RATE*nH);
}

// run a step of the derivatives
void derivs(double x, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i += 2) {		
		double *sub_vec_nHx = vec_nHx + (i-1);
		
		vec_dnHxdt[i] = dnH2(sub_vec_nHx);
		vec_dnHxdt[i+1] = dnH_p(sub_vec_nHx);
	}
}

double (*derivs_arr[])(double[]) = { &dnH2, &dnH_p };

void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar)
{
	for (int i = 1; i <= nvar; i += 2) {
		dfdx[i] = 0.0;
		dfdx[i+1] = 0.0;
		
		double *sub_vec_nHx = vec_nHx + (i-1);
		jacobian(derivs_arr, 2, sub_vec_nHx, 1e-6, dfdy, nvar);
	}
}

double k1(double temp, double temp_grain) { 
	double fa = 1/(1 + 1e4*exp(-600/temp_grain));
	double t2 = temp / 100;
	double t_grain2 = temp_grain / 100;
	
	return STEP_TIME * 3e-17 * ((sqrt(t2) * fa) / (1+0.4*sqrt(t2+t_grain2) + 0.2*t2 + 0.08*t2*t2));
}
double k2(double temp, double nH2, double nH_p) {
	double kH = 3.52e-9*exp(-4.39e4/temp);
	double kL = 6.67e-12*sqrt(temp)*exp(-(1+63590/temp)); // Glover et al.
	
	if (kL < 1e-16) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(10, 3.0 - 0.416*log_temp - 0.327*log_temp*log_temp);
	double ncr_H2 = pow(10, 4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	double n = N_tot * ((getxH(nH2,nH_p)/ncr_H) + (getxH2(nH2)/ncr_H2));
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double k3(double temp, double nH2, double nH_p) { 
	double kH = 1.3e-9*exp(-5.33e4/temp);	
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	if (kL == 0) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(10, 3.0 - 0.416*log_temp - 0.327*log_temp*log_temp);
	double ncr_H2 = pow(10, 4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	double n = N_tot * ((getxH(nH2,nH_p)/ncr_H) + (getxH2(nH2)/ncr_H2));
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double k6(double temp) {
	// Drains Physics of the Interstellar and intergalactic medium
	return STEP_TIME * 5.466e-9*1.07*sqrt(temp/1.0e4)*exp(-13.6/(8.6173e-5*temp));
}
double k7(double temp) {
	return STEP_TIME * 2.54e-13*pow(temp/1.0e4, -0.8163 - 0.0208*log(temp/1.0e4));
}
double k8(double temp, double nH_p) {
	double G0 = 1.13,
		   ne = getne(nH_p),
		   psi = (G0*sqrt(temp)) / ne;
	return STEP_TIME * (1e-14*12.25)/(1 + 8.074e-6*pow(psi, 1.378)*(1 + 5.087e2*pow(temp,1.586e-2)*pow(psi,-0.4723 - 1.102e-5*log(temp))));
}