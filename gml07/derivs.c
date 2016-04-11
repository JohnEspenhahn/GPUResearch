#include "derivs.h"

/* private */ int TEMP = TEMP_INIT, GRAIN_TEMP = GRAIN_TEMP_INIT;

double getTemp() { return TEMP; }
void setTemp(int t) { if (t > 0) TEMP = t; }

double getGrainTemp() { return GRAIN_TEMP; }
void setGrainTemp(int t) { if (t > 0) GRAIN_TEMP = t; }

/** Takes mass densities and converts to number density */
double getnH(double pH2, double pH_p) { return N_H_tot - 2*getnH2(pH2) - getnH_p(pH_p); }

/** Takes mass density of H2 and converts to number density */
double getnH2(double pH2) { return pH2 / (mu_h2*M_h); }

/** Takes mass density of H+ and converts to number density */
double getnH_p(double pH_p) { return pH_p / (mu_h*M_h); }

/** Takes mass density of H2 and converts to abundance */
double getxH2(double pH2) { return (2*getnH2(pH2)) / N_tot; }

/** Takes mass density of H_p and converts to abundance */
double getxH_p(double pH_p) { return getnH_p(pH_p) / N_tot; }

/** Derives the abundance of atomic hydrogen based on mass densities of H2 and H+ */
double getxH(double pH2, double pH_p) { 
	return xH_tot - getxH2(pH2) - getxH_p(pH_p); 
}

/** Derives the mass density of atomic hydrogen based on mass densities of H2 and H+ */
double getpH(double pH2, double pH_p) { return getnH(pH2, pH_p)*(mu_h*M_h); };

/** Derives the abundance of e- based on the mass densities of H2 and H+ */
double getxe(double pH_p) { 
	return getxH_p(pH_p) + xC + xSi; 
}

/** Derives the number density of e- based on the mass density of H+ */
double getne(double pH_p) {
	return getxe(pH_p) * N_tot;
}

double dpH2(double vec_pHx[]) {
	double pH2 = vec_pHx[1]
		 , pH_p = vec_pHx[2]
		 , nH = getnH(pH2, pH_p)
		 , nH2 = getnH2(pH2)
		 , nH_p = getnH_p(pH_p);
	
	return M_h*mu_h2 * (k1(TEMP, GRAIN_TEMP)*nH*nH - k2(TEMP,nH2,nH_p)*nH2*nH - k3(TEMP,nH2,nH_p)*nH2*nH2 - SELF_SHIELDING*nH2);
}

double dpH_p(double vec_pHx[]) {
	double pH2  = vec_pHx[1]
	     , pH_p = vec_pHx[2]
		 , nH   = getnH(pH2, pH_p)
		 , nH_p = getnH_p(pH_p)
		 , ne   = getne(pH_p);
	
	return M_h*mu_h * (k6(TEMP)*nH*ne - k7(TEMP)*nH_p*ne - k8(TEMP,ne)*nH_p*ne + COSMIC_RAY_RATE*nH);
}

// run a step of the derivatives
void derivs(double x, int nvar, double vec_pHx[], double vec_dpHxdt[]) {
	for (int i = 1; i <= nvar; i += 2) {		
		double *sub_vec_pHx = vec_pHx + (i-1);
		
		// pH2
		vec_dpHxdt[i] = dpH2(sub_vec_pHx);
		
		// pH_p
		vec_dpHxdt[i+1] = dpH_p(sub_vec_pHx);
	}
}

double (*derivs_arr[])(double[]) = { &dpH2, &dpH_p };

void jacobn(double x, double vec_pHx[], double dfdx[], double **dfdy, int nvar)
{
	for (int i = 1; i <= nvar; i += 2) {
		dfdx[i] = 0.0;
		dfdx[i+1] = 0.0;
		
		/*
		double *sub_vec_pHx = vec_pHx + (i-1);
		jacobian(derivs_arr, 2, sub_vec_pHx, 1e-32, dfdy, nvar);
		*/
		
		double pH2 = vec_pHx[i]
			 , pH_p = vec_pHx[i+1]
			 , pH = getpH(pH2, pH_p)
			 , nH = getnH(pH2, pH_p);
		
		// dpH2 / dpH2
		dfdy[i][1] = k1(TEMP, GRAIN_TEMP)*pH*(-1/(mu_h*M_h)) 
						- k2(TEMP, pH2, pH_p)*nH 
						- k2(TEMP, pH2, pH_p)*pH2*(-2/(mu_h2*M_h)) 
						- (4*k3(TEMP, pH2, pH_p)*pH2)/(mu_h2*M_h)
						+ SELF_SHIELDING;
						
		// dpH2 / dpH+
		dfdy[i][2] = -2*k1(TEMP, GRAIN_TEMP)*nH 
						- (k2(TEMP, pH2, pH_p)*pH2)/(mu_h*M_h);
		
		// dpH+ / dpH2
		dfdy[i+1][1] = -k6(TEMP)*(2*mu_h/mu_h2)*(pH_p/(mu_h*M_h) + nSi + nC)
						- COSMIC_RAY_RATE*(2*mu_h/mu_h2);
		
		// dpH+ / dpH+	
		dfdy[i+1][2] = k6(TEMP)*(N_H_tot - nC - nSi - 2*pH2/(mu_h2*M_h) - 2*pH_p/(mu_h*M_h)) 
						- (k7(TEMP) + k8(TEMP,pH_p))*(2*pH_p/(mu_h*M_h) + nSi + nC)
						- COSMIC_RAY_RATE;
	}
}

double k1(double temp, double temp_grain) { 
	double fa = 1/(1 + 1e4*exp(-600/temp_grain));
	double t2 = temp / 100;
	double t_grain2 = temp_grain / 100;
	
	return STEP_TIME * 3e-17 * ((sqrt(t2) * fa) / (1+0.4*sqrt(t2+t_grain2) + 0.2*t2 + 0.08*t2*t2));
}
double k2(double temp, double pH2, double pH_p) {
	double kH = 3.52e-9*exp(-4.39e4/temp);
	double kL = 6.67e-12*sqrt(temp)*exp(-(1+63590/temp)); // Glover et al.
	
	if (kL < 1e-16) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(10, 3.0 - 0.416*log_temp - 0.327*log_temp*log_temp);
	double ncr_H2 = pow(10, 4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	double n = N_tot * ((getxH(pH2,pH_p)/ncr_H) + (getxH2(pH2)/ncr_H2));
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double k3(double temp, double pH2, double pH_p) { 
	double kH = 1.3e-9*exp(-5.33e4/temp);	
	double kL = 1.18e-10*exp(-6.95e4/temp);
	
	if (kL == 0) return 0;
	
	double log_temp = log10(temp/1.0e4);
	double ncr_H = pow(10, 3.0 - 0.416*log_temp - 0.327*log_temp*log_temp);
	double ncr_H2 = pow(10, 4.845 - 1.3*log_temp + 1.62*log_temp*log_temp);
	
	double n = N_tot * ((getxH(pH2,pH_p)/ncr_H) + (getxH2(pH2)/ncr_H2));
	return STEP_TIME * pow(kH, n/(1+n)) * pow(kL, 1/(1+n));
}
double k6(double temp) {
	// Drains Physics of the Interstellar and intergalactic medium
	return STEP_TIME * 5.466e-9*1.07*sqrt(temp/1.0e4)*exp(-13.6/(8.6173e-5*temp));
}
double k7(double temp) {
	return STEP_TIME * 2.54e-13*pow(temp/1.0e4, -0.8163 - 0.0208*log(temp/1.0e4));
}
double k8(double temp, double pH_p) {
	double G0 = 1.13,
		   ne = getne(pH_p),
		   psi = (G0*sqrt(temp)) / ne;
	return STEP_TIME * (1e-14*12.25)/(1 + 8.074e-6*pow(psi, 1.378)*(1 + 5.087e2*pow(temp,1.586e-2)*pow(psi,-0.4723 - 1.102e-5*log(temp))));
}