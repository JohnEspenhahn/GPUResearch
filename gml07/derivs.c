#include "derivs.h"

/* private */ int TEMP = TEMP_INIT, GRAIN_TEMP = GRAIN_TEMP_INIT;

double getTemp() { return TEMP; }
void setTemp(int t) { if (t > 0) TEMP = t; }

double getGrainTemp() { return GRAIN_TEMP; }
void setGrainTemp(int t) { if (t > 0) GRAIN_TEMP = t; }

/** Takes mass density of H2 and converts to number density */
double getnH2(double pH2) { return pH2 / (2.36 * M_h); }

/** Takes mass density of H+ and converts to number density */
double getnH_p(double pH_p) { return pH_p / (1.36 * M_h); }

/** Takes mass density of H2 and converts to abundance */
double getxH2(double pH2) { return (2*getnH2(pH2)) / N_TOT; }

/** Takes mass density of H_p and converts to abundance */
double getxH_p(double pH_p) { return getnH_p(pH_p) / N_TOT; }

/** Gets the total adundance of hydrogen nuclei */
double getxH_tot() { return xH_tot; }

/** Derives the mass density of atomic hydrogen based on mass densities of H2 and H+ */
double getpH(double pH2, double pH_p) { return (getxH(pH2, pH_p) * N_TOT) * 1.36 * M_h; };

/** Derives the abundance of atomic hydrogen based on mass densities of H2 and H+ */
double getxH(double pH2, double pH_p) { 
	return xH_tot - getxH2(pH2) - getxH_p(pH_p); 
}

/** Derives the abundance of e- based on the mass densities of H2 and H+ */
double getxe(double pH2, double pH_p) { 
	return getxH_p(pH_p) + xC + xSi; 
}

/** Derives the number density of e- based on the mass densities of H2 and H+ */
double getne(double pH2, double pH_p) {
	return getxe(pH2, pH_p) * N_TOT;
}

// run a step of the derivatives
void derivs(double t, int nvar, double vec_pHx[], double vec_dpHxdt[]) {
	for (int i = 1; i <= nvar; i += 2) {
		double pH2 = vec_pHx[i];
		double pH_p = vec_pHx[i+1];
		double pH = getpH(pH2, pH_p); 
		double ne = getne(pH2, pH_p);
		
		// pH2
		vec_dpHxdt[i] = k1(TEMP, GRAIN_TEMP)*pH*pH - k2(TEMP,pH2,pH_p)*pH2*pH - k3(TEMP,pH2,pH_p)*pH2*pH2 - SELF_SHIELDING*pH2;
		
		// pH_p
		vec_dpHxdt[i+1] = k6(TEMP)*pH*ne - k7(TEMP)*pH_p*ne - k8(TEMP,pH2,pH_p)*pH_p*ne + COSMIC_RAY_RATE*pH;
	}
}

void jacobn(double x, double vec_pHx[], double dfdx[], double **dfdy, int nvar)
{
	for (int i = 1; i <= nvar; i += 2) {
		dfdx[i] = 0.0;
		dfdx[i+1] = 0.0;
		
		double pH2 = vec_pHx[i],
			   xH2 = getxH2(pH2),
			   pH_p = vec_pHx[i+1],
			   pH = getpH(pH2, pH_p),
			   ne = getne(pH2, pH_p),
			   xH = getxH(pH2,pH_p);
		
		// pH2
		dfdy[i][1] = -k2(TEMP,xH2,xH)*pH - 2*k3(TEMP,pH2,pH_p)*pH2 - dk3ndH2(TEMP,pH2)*pH2*pH2 - SELF_SHIELDING; // need dk3ndH2 if k3 is actually function of pH2
		dfdy[i][2] = 0;
		
		// pH_p
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
double k2(double temp, double xH2, double xH) {
	double kH = 3.52e-9*exp(-4.39e4/temp);
	double kL = 1.12e-10*exp(-7.035e4/temp);
	
	if (kL == 0) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(3.0 - 0.416*log_temp - 0.327*log_temp*log_temp, 10);
	double ncr_H2 = pow(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp, 10);
	
	double n = N_TOT * ((xH/ncr_H) + (xH2/ncr_H2));
	
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double k3(double temp, double pH2, double pH_p) { 
	double kH = 1.3e-9*exp(-5.33e4/temp);	
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	if (kL == 0) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(3.0 - 0.416*log_temp - 0.327*log_temp*log_temp, 10);
	double ncr_H2 = pow(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp, 10);
	
	double n = N_TOT * ((getxH(pH2,pH_p)/ncr_H) + (getxH2(pH2)/ncr_H2));
	
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double dk3ndH2(double temp, double pH2) {
	double kH = 1.3e-9*exp(-5.33e4/temp);	
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	if (kL == 0) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr = exp(4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	return STEP_TIME * k3(temp, pH2) * ((ncr*log10(kH/kL)) / ((ncr+pH2)*(ncr+pH2)));
}
double k6(double temp) {
	return STEP_TIME * 5.466e-9*1.07*sqrt(temp/1.0e4)*exp(-13.6/(8.6173e-5*temp));
}
double k7(double temp) {
	return STEP_TIME * 2.54e-13*pow(temp/1.0e4, -0.8163 - 0.0208*log(temp/1.0e4));
}
double k8(double temp, double pH2, double pH_p) {
	double G0 = 1.13;
	double psi = (G0*sqrt(temp)) / getne(pH2, pH_p);
	return STEP_TIME * (1e-14*12.25)/(1 + 8.074e-6*pow(psi, 1.378)*(1 + 5.087e2*pow(temp,1.586e-2)*pow(psi,-0.4723 - 1.102e-5*log(temp))));
}