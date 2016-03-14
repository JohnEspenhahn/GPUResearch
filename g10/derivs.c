#include "derivs.h"

int T = TEMP_INIT, GRAIN_TEMP = GRAIN_TEMP_INIT;

// run a step of the derivatives
void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i += 2) {
		double nH2 = vec_nHx[i];
		double nH_p = vec_nHx[i+1];
		double nH = getnH(nH2, nH_p); 
		double ne = getne(nH2, nH_p);
		
		// pH_p
		vec_dnHxdt[i] = k4(t)*pH*pH2_p + k11(t)*pe*pH + k18(t)*pH*pHe_p + k24(t)*pO_p*pH + k28(t)*pC_p*pH + k89(t)*pH2*pHe_p + k97(t)*pH2O*pHe_p + k106(t)*pCO_p*pH - k3(t)*pH*pH_p - k5(t)*pH_m*pH_p - k7(t)*pH_p*pH2 - k12(t)*pe*pH_p - k15(t)*pH_m*pH_p - k19(t)*pH_p*pHe - k25(t)*pH_p*pO - k27(t)*pC*pH_p - k62(t)*pH_p*pCH2 - k90(t)*pCH*pH_p - k91(t)*pH_p*pCH2 - k94(t)*pH_p*pOH - k96(t)*pH2O*pH_p - k100(t)*pO2*pH_p - k107(t)*pC_m*pH_p - k108(t)*pH_p*pO_m - k141(t)*pH_p*pH2;
		
		// pH2
		vec_dnHxdt[i+1] = k2(t)*pH_m*pH + k4(t)*pH*pH2_p + k35(t)*pCH*pH + k39(t)*pH*pCH2 + k41(t)*pCH2*pO + k44(t)*pH*pOH + k49(t)*pH2O*pH + k55(t)*pH*pH3_p + k57(t)*pC*pH3_p + k59(t)*pH*pCH_p + k62(t)*pH_p*pCH2 + k63(t)*pH*pCH2_p + k66(t)*pH*pCH3_p + k67(t)*pCH3_p*pO + k71(t)*pH3_p*pO + k72(t)*pOH*pH3_p + k76(t)*pH2O*pH3_p + k79(t)*pC*pH3O_p + k84(t)*pCO*pH3_p + k85(t)*pCO*pH3_p + k92(t)*pCH2*pHe_p + k110(t)*pe*pH3_p + k115(t)*pe*pCH2_p + k117(t)*pe*pCH3_p + k121(t)*pe*pH2O_p + k124(t)*pe*pH3O_p + k126(t)*pe*pH3O_p + k154(t)*3*pH + k155(t)*2*pH*pH2 + k156(t)*2*pH*pHe + k165(t)*pH*pH - k7(t)*pH_p*pH2 - k8(t)*pe*pH2 - k9(t)*pH*pH2 - k10(t)*pH2*pH2 - k30(t)*pH2*pHe - k34(t)*pC*pH2 - k36(t)*pCH*pH2 - k43(t)*pH2*pO - k45(t)*pOH*pH2 - k51(t)*pO2*pH2 - k54(t)*pH2_p*pH2 - k58(t)*pC_p*pH2 - k60(t)*pCH_p*pH2 - k64(t)*pH2*pCH2_p - k69(t)*pO_p*pH2 - k74(t)*pH2*pOH_p - k75(t)*pH2O_p*pH2 - k88(t)*pH2*pHe_p - k89(t)*pH2*pHe_p - k136(t)*pC_m*pH2 - k139(t)*pH2*pO_m - k141(t)*pH_p*pH2 - k144(t)*pC*pH2 - k148(t)*pC_p*pH2;
		
		// pHe_p
		vec_dnHxdt[i+2] = k16(t)*pe*pHe + k19(t)*pH_p*pHe - k17(t)*pe*pHe_p - k18(t)*pH*pHe_p - k26(t)*pHe_p*pO - k29(t)*pC*pHe_p - k88(t)*pH2*pHe_p - k89(t)*pH2*pHe_p - k92(t)*pCH2*pHe_p - k93(t)*pHe_p*pC2 - k95(t)*pOH*pHe_p - k97(t)*pH2O*pHe_p - k98(t)*pH2O*pHe_p - k99(t)*pH2O*pHe_p - k101(t)*pO2*pHe_p - k102(t)*pO2*pHe_p - k104(t)*pCO*pHe_p - k105(t)*pCO*pHe_p - k109(t)*pH_m*pHe_p;
		
		// pC_p
		vec_dnHxdt[i+3] = k22(t)*pC*pe + k27(t)*pC*pH_p + k29(t)*pC*pHe_p + k59(t)*pH*pCH_p + k92(t)*pCH2*pHe_p + k93(t)*pHe_p*pC2 + k103(t)*pO2_p*pC + k104(t)*pCO*pHe_p - k20(t)*pC_p*pe - k28(t)*pC_p*pH - k58(t)*pC_p*pH2 - k73(t)*pC_p*pOH - k77(t)*pH2O*pC_p - k78(t)*pH2O*pC_p - k80(t)*pO2*pC_p - k81(t)*pO2*pC_p - k147(t)*pC_p*pH - k148(t)*pC_p*pH2 - k149(t)*pC_p*pO - k159(t)*pC_p*pM*pO;
		
		// pO_p
		vec_dnHxdt[i+4] = k23(t)*pe*pO + k25(t)*pH_p*pO + k26(t)*pHe_p*pO + k81(t)*pO2*pC_p + k95(t)*pOH*pHe_p + k102(t)*pO2*pHe_p + k105(t)*pCO*pHe_p - k21(t)*pO_p*pe - k24(t)*pO_p*pH - k68(t)*pO_p*pC2 - k69(t)*pO_p*pH2 - k160(t)*pO_p*pC*pM;
		
		// pOH
		vec_dnHxdt[i+5] = k43(t)*pH2*pO + k49(t)*pH2O*pH + k50(t)*pO2*pH + k51(t)*pO2*pH2 + k53(t)*pH*pCO + k82(t)*pO2*pCH2_p + k97(t)*pH2O*pHe_p + k122(t)*pe*pH2O_p + k124(t)*pe*pH3O_p + k125(t)*pe*pH3O_p + k130(t)*pHCO_p*pe + k133(t)*pH_m*pO + k138(t)*pH*pO_m + k151(t)*pH*pO + k161(t)*pH*pM*pO - k31(t)*pH*pOH - k44(t)*pH*pOH - k45(t)*pOH*pH2 - k46(t)*pC*pOH - k47(t)*pOH*pO - k48(t)*2*pOH - k72(t)*pOH*pH3_p - k73(t)*pC_p*pOH - k94(t)*pH_p*pOH - k95(t)*pOH*pHe_p - k134(t)*pH_m*pOH - k153(t)*pH*pOH - k162(t)*pH*pOH*pM;
		
		// pH2O
		vec_dnHxdt[i+6] = k45(t)*pOH*pH2 + k48(t)*pOH*pOH + k123(t)*pe*pH3O_p + k134(t)*pH_m*pOH + k139(t)*pH2*pO_m + k153(t)*pH*pOH + k162(t)*pH*pOH*pM - k49(t)*pH2O*pH - k76(t)*pH2O*pH3_p - k77(t)*pH2O*pC_p - k78(t)*pH2O*pC_p - k87(t)*pHCO_p*pH2O - k96(t)*pH2O*pH_p - k97(t)*pH2O*pHe_p - k98(t)*pH2O*pHe_p - k99(t)*pH2O*pHe_p;
		
		//pCO
		vec_dnHxdt[i+7] = k38(t)*pCH*pO + k40(t)*pCH2*pO + k41(t)*pCH2*pO + k42(t)*pC2*pO + k46(t)*pC*pOH + k52(t)*pO2*pC + k81(t)*pO2*pC_p + k86(t)*pHCO_p*pC + k87(t)*pHCO_p*pH2O + k106(t)*pCO_p*pH + k129(t)*pHCO_p*pe + k131(t)*pe*pHOC_p + k137(t)*pC_m*pO + k140(t)*pC*pO_m + k146(t)*pC*pO + k158(t)*pC*pM*pO - k53(t)*pH*pCO - k84(t)*pCO*pH3_p - k85(t)*pCO*pH3_p - k104(t)*pCO*pHe_p - k105(t)*pCO*pHe_p;
		
		// pC2
		vec_dnHxdt[i+8] = k37(t)*pC*pCH + k145(t)*2*pC + k157(t)*2*pC*pM - k42(t)*pC2*pO - k68(t)*pO_p*pC2 - k93(t)*pHe_p*pC2;
		
		// pO2
		vec_dnHxdt[i+9] = k47(t)*pOH*pO + k103(t)*pO2_p*pC + k152(t)*2*pO + k163(t)*pM*2*pO - k50(t)*pO2*pH - k51(t)*pO2*pH2 - k52(t)*pO2*pC - k80(t)*pO2*pC_p - k81(t)*pO2*pC_p - k82(t)*pO2*pCH2_p - k100(t)*pO2*pH_p - k101(t)*pO2*pHe_p - k102(t)*pO2*pHe_p;
		
		// pHCO_p
		vec_dnHxdt[i+10] = k32(t)*pHOC_p*pH2 + k33(t)*pHOC_p*pCO + k65(t)*pCH2_p*pO + k67(t)*pCH3_p*pO + k77(t)*pH2O*pC_p + k78(t)*pH2O*pC_p + k79(t)*pC*pH3O_p + k82(t)*pO2*pCH2_p + k85(t)*pCO*pH3_p + k164(t)*pCH*pO - k86(t)*pHCO_p*pC - k87(t)*pHCO_p*pH2O - k129(t)*pHCO_p*pe - k130(t)*pHCO_p*pe;
		
		// pCH
		vec_dnHxdt[i+11] = k34(t)*pC*pH2 + k39(t)*pH*pCH2 + k113(t)*pe*pCH2_p + k117(t)*pe*pCH3_p + k118(t)*pe*pCH3_p + k132(t)*pC*pH_m + k135(t)*pC_m*pH + k143(t)*pC*pH - k35(t)*pCH*pH - k36(t)*pCH*pH2 - k37(t)*pC*pCH - k38(t)*pCH*pO - k90(t)*pCH*pH_p - k164(t)*pCH*pO;
		
		// pCH2
		vec_dnHxdt[i+12] = k36(t)*pCH*pH2 + k116(t)*pe*pCH3_p + k136(t)*pC_m*pH2 + k144(t)*pC*pH2 - k39(t)*pH*pCH2 - k40(t)*pCH2*pO - k41(t)*pCH2*pO - k62(t)*pH_p*pCH2 - k91(t)*pH_p*pCH2 - k92(t)*pCH2*pHe_p;
		
		// pCH3_p
		vec_dnHxdt[i+13] = k64(t)*pH2*pCH2_p - k66(t)*pH*pCH3_p - k67(t)*pCH3_p*pO - k116(t)*pe*pCH3_p - k117(t)*pe*pCH3_p - k118(t)*pe*pCH3_p;
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
		dfdy[i][1] = -k2(T,nH)*nH - 2*k3(T,nH)*nH2 - dk3ndH2(T,nH2)*nH2*nH2 - SELF_SHIELDING; // need dk3ndH2 if k3 is actually function of nH2
		dfdy[i][2] = 0;
		
		// nH_p
		dfdy[i+1][1] = 0;
		dfdy[i+1][2] = -k7(T)*ne - k8(T)*ne;
	}
}

double k1(double t) {
	double lt = log10(t);
	if (t <= 6000) return pow(-17.845 + 0.762*lt + 0.1523*lt*lt - 0.03274*lt*lt*lt, 10);
	else return pow(-16.420 + 0.1998*lt*lt - 5.447e-3*lt*lt*lt*lt + 4.0415e-5*lt*lt*lt*lt*lt*lt, 10);
}
double k2(double t) {
	if (t <= 300) return 1.5e-9;
	else return 4.0e-9*pow(t, -0.17);
}
double k3(double t) {
	double lt = log10(t);
	return pow(-19.38 - 1.523*lt + 1.118*lt*lt - 0.1269*lt*lt*lt, 10);
}
double k4(double t) {
	return 6.4e-10;
}
double k5(double t) {
	return 2.4e-6*pow(t, -0.5)*(1.0 + t/20000.0);
}
double k6(double t) {
	if (t <= 617) return 1.0e-8;
	else return 1.32e-6*pow(t, -0.76);
}
double k7(double t) {
	double lt = log(t); // ln
	return exp(-21237.15/t)*(-3.3232183e-7 + 3.3735382e-7*lt - 1.4491368e-7*lt*lt + 3.4172805e-8*lt*lt*lt - 4.7813720E-9*lt*lt*lt*lt + 3.9731542e-10*lt*lt*lt*lt*lt - 1.8171411e-11*ln*ln*ln*ln*ln*ln + 3.5311932e-13*ln*ln*ln*ln*ln*ln*ln);
}
double k8(double t) {
	return 3.73e-9*pow(t,0.1121)*exp(-99430.0/t);
}
double k9(double t, double nH) {
	double kl = 6.67e-12*pow(t)*exp(-(1+63590/t)),
		   kh = 3.52e-9*exp(-43900/t),
		   logt = log10(t/10000),
		   ncr = nH / pow(3.0 - 0.416*logt - 0.327*logt*logt, 10);
	return kh*pow(kh/kl,-1/(ncr+1));
}
double k10(double t, double nH2) {
	double kl = (5.996e-30*pow(t,4.1881)*exp(-54657.4/t))/pow(1.0 + 6.761e-6 * t, 5.6881),
		   kh = 1.3e-9*exp(-53300/t),
		   logt = log10(t/10000),
		   ncr = nH2 / pow(4.845 - 1.3*logt + 1.62*logt*logt, 10);
	return kh*pow(kh/kl,-1/(ncr+1));
}
double k11(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-32.71396786 + 13.5365560*lt - 5.73932875*lt*lt + 1.56315498*lt*lt*lt - 2.87705600e−1*lt*lt*lt*lt + 3.48255977e−2*lt*lt*lt*lt*lt - 2.63197617e−3*lt*lt*lt*lt*lt*ln + 1.11954395e−4*lt*lt*lt*lt*lt*ln*ln - 2.03914985e−6*lt*lt*lt*lt*lt*ln*ln*ln);
}
double k12(double t) {
	return 1.269e-13*pow(315614/t, 1.503)*pow(1+pow(6046255/t,0.470),-1.923);
}
double k13(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-1.801849334e1 + 2.36085220*lt - 2.82744300e−1*lt*lt + 1.62331664e−2*lt*lt*lt - 3.36501203e−2*lt*lt*lt*lt + 1.17832978e−2*lt*lt*lt*lt*lt - 1.65619470e−3*lt*lt*lt*lt*lt*lt + 1.06827520e−4*lt*lt*lt*lt*lt*lt*lt - 2.63128581e−6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k14(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	if (te <= 0.1) return 2.5634e-9*pow(te,1.78186);
	else return exp(-2.0372609e1 + 1.13944933e0*lt - 1.4210135e−1*lt*lt + 8.4644554e−3*lt*lt*lt - 1.4327641e−3*lt*lt*lt*lt + 2.0122503e−4*lt*lt*lt*lt*lt + 8.6639632e−5*lt*lt*lt*lt*lt*lt - 2.5850097e−5*lt*lt*lt*lt*lt*lt*lt + 2.4555012e−6*lt*lt*lt*lt*lt*lt*lt*lt - 8.0683825e−8*lt*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k15(double t) {
	if (t <= 8000) return 6.9e-9*pow(t,-0.35);
	else return 9.6e-7*pow(t,-0.90);
}
double k16(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-4.409864886e1 + 2.391596563e1*lt - 1.07532302e1*lt*lt + 3.05803875e0*lt*lt*lt - 5.6851189e−1*lt*lt*lt*lt + 6.79539123e−2*lt*lt*lt*lt*lt - 5.0090561e−3*lt*lt*lt*lt*lt*lt + 2.06723616e−4*lt*lt*lt*lt*lt*lt*lt - 3.64916141e−6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k17(double t) {
	return 
}
double k18(double t) {
	return 1.25e-15*pow(t/300,0.25);
}
double k19(double t) {
	if (t <= 10000) return 1.26e-9*pow(t,-0.75)*exp(-127500/t);
	else return 4.0e-37*pow(t,4.74);
}
double k20(double t) {
	if (t <= 7950) return 4.67e-12*pow(t/300,-0.6);
	else if (t <= 21140) return 1.23e-17*pow(t/300,2.49)*exp(21845.6/t);
	else return 9.62e-8*pow(t/300,-1.37)*exp(-115786.2/t);
}
double k21(double t) {
	if (t <= 400) return 1.30e-10*pow(t,-0.64);
	else return 1.41e-10*pow(t,-0.66) + 7.4e-4*pow(t,-1.5)*exp(-175000/t)*(1.0 + 0.062*exp(-145000/t));
}
double k22(double t) {
	double te = t * 8.621738e-5;
	double u = 11.26/te;
	return (6.85e-8/(0.193 + u)) * pow(u,0.25) * exp(-u);
}
double k23(double t) {
	double te = t * 8.621738e-5;
	double u = 13.6/te;
	return (3.59e-8/(0.073 + u)) * pow(u,0.34) * exp(-u);
}
double k24(double t) {
	return 4.99e-11 * pow(t,0.405) + 7.54e-10*pow(t,-0.458);
}
double k25(double t) {
	return (1.08e−11*pow(t,0.517) + 4.00e−10*pow(t,0.00669)) * exp(-227/t);
}
double k26(double t) {
	return 4.991e-15*pow(t/10000,0.3794)*exp(-t/1121000) + 2.78e-15*pow(t/10000,-0.2163)*exp(t/815800);
}
double k27(double t) {
	return 3.9e-16*pow(t,0.213);
}
double k28(double t) {
	return 6.08e-14*pow(t/10000,1.96)*exp(-170000/t);
}
double k29(double t) {
	if (t <= 200) return 8.58e-17*pow(t,0.757);
	else if (t <= 2000) return 3.25e-17*pow(t,0.968);
	else return 2.77e-19*pow(t,1.597);
}
double k30(double t, double nH2) {
	double logt = log10(t),
		   kl = pow(-27.029 + 3.801*logt - 29487/t, 10),
		   kh = pow(-2.729 - 1.75*logt - 23474/t, 10),
		   ncr = pow(nH2 / pow(5.0792*(1.0 - 1.23e-5*(t-2000)), 10), 1.07);
	return kh*pow(kh/kl,-1/(ncr+1));
}
double k31(double t) {
	return 6.0e-9*exp(-50900/t);
}
double k32(double t) {
	return 3.8e-10;
}
double k33(double t) {
	return 4.0e-10;
}
double k34(double t) {
	return 6.64e-10*exp(-11700/t);
}
double k35(double t) {
	return 1.31e-10*exp(-80/t);
}
double k36(double t) {
	return 5.46e-10*exp(-1943/t);
}
double k37(double t) {
	return 6.59e-11;
}
double k38(double t) {
	if (t <= 2000) return 6.6e-11;
	else if (t > 2000) return 1.02e-10*exp(-914/t);
}
double k39(double t) {
	return 6.64e-11;
}
double k40(double t) {
	return 1.33e-10;
}
double k41(double t) {
	return 8.0e-11;
}
double k42(double t) {
	if (t <= 300) return 5.0e-11*sqrt(t/500);
	else return 5.0e-11*pow(t/300,0.757);
}
double k43(double t) {
	return 3.14e-13*pow(t/300,2.7)*exp(-3150/t);
}
double k44(double t) {
	return 6.99e-14*pow(t/300,2.8)*exp(-1950/t);
}
double k45(double t) {
	return 2.05e-12*pow(t/300,1.52)*exp(-1736/t);
}
double k46(double t) {
	return 1e-10;
}
double k47(double t) {
	if (t <= 261) return 3.5e-11;
	else 1.77e-11*exp(178/t);
}
double k48(double t) {
	return 1.65e-12*pow(t/300,1.14)*exp(-50/t);
}
double k49(double t) {
	return 1.59e-11*pow(t/300,1.2)*exp(-9610/t);
}
double k50(double t) {
	return 2.61e-10*exp(-8156/t);
}
double k51(double t) {
	return 3.16e-10*exp(-21890/t);
}
double k52(double t) {
	if (t <= 295) return 4.7e-11*pow(t/300,-0.34);
	else return 2.48e-12*pow(t/300,1.54)*exp(613/t);
}
double k53(double t) {
	return 1.1e-10*sqrt(t/300)*exp(-77700/t);
}
double k54(double t) {
	return 2.24e-9*pow(t/300,0.042)*exp(-t/466000);
}
double k55(double t) {
	return 7.7e-9*exp(-17560/t);
}
double k56(double t) {
	return 2.4e-9;
}
double k57(double t) {
	return 2.0e-9;
}
double k58(double t) {
	return 1e-10*exp(-4640/t);
}
double k59(double t) {
	return 7.5e-10;
}
double k60(double t) {
	return 1.2e-9;
}
double k61(double t) {
	return 3.5e-10;
}
double k62(double t) {
	return 1.4e-9;
}
double k63(double t) {
	return 1e-9*exp(-7080/t);
}
double k64(double t) {
	return 1.6e-9;
}
double k65(double t) {
	return 7.5e-10;
}
double k66(double t) {
	return 7.0e-10*exp(-10560/t);
}
double k67(double t) {
	return 4.0e-10;
}
double k68(double t) {
	return 4.8e-10;
}
double k69(double t) {
	return 1.7e-9;
}
double k70(double t) {
	return 1.5e-9;
}
double k71(double t) {
	return 8.4e-10;
}
double k72(double t) {
	return 1.3e-9;
}
double k73(double t) {
	return 7.7e-10;
}
double k74(double t) {
	return 1.01e-9;
}
double k75(double t) {
	return 6.4e-10;
}
double k76(double t) {
	return 5.9e-9;
}
double k77(double t) {
	return 9.0e-10;
}
double k78(double t) {
	return 1.8e-9;
}
double k79(double t) {
	return 1e-11;
}
double k80(double t) {
	return 3.8e-10;
}
double k81(double t) {
	return 6.2e-10;
}
double k82(double t) {
	return 9.1e-10;
}
double k83(double t) {
	return 5.2e-11;
}
double k84(double t) {
	return 2.7e-11;
}
double k85(double t) {
	return 1.7e-9;
}
double k86(double t) {
	return 1.1e-9;
}
double k87(double t) {
	return 2.5e-9;
}
double k88(double t) {
	return 7.2e-15;
}
double k89(double t) {
	return 3.7e-14*exp(-35/t);
}
double k90(double t) {
	return 1.9e-9;
}
double k91(double t) {
	return 1.4e-9;
}
double k92(double t) {
	return 7.5e-10;
}
double k93(double t) {
	return 1.6e-9;
}
double k94(double t) {
	return 2.1e-9;
}
double k95(double t) {
	return 1.1e-9;
}
double k96(double t) {
	return 6.9e-9;
}
double k97(double t) {
	return 2.04e-10;
}
double k98(double t) {
	return 2.86e-10;
}
double k99(double t) {
	return 6.05e-11;
}
double k100(double t) {
	return 2e-9;
}
double k101(double t) {
	return 3.3e-11;
}
double k102(double t) {
	return 1.1e-9;
}
double k103(double t) {
	return 5.2e-11;
}
double k104(double t) {
	return 1.4e-9*pow(t/300,-0.5);
}
double k105(double t) {
	return 1.4e-16*pow(t/300,-0.5);
}
double k106(double t) {
	return 7.5e-10;
}
double k107(double t) {
	return 2.3e-7*pow(t/300,-0.5);
}
double k108(double t) {
	return 2.3e-7*pow(t/300,-0.5);
}
double k109(double t) {
	return 2.32e-7*pow(t/300,-0.52)*exp(t/22400);
}
double k110(double t) {
	return 2.34e-8*pow(t/300,-0.52);
}
double k111(double t) {
	return 4.36e-8*pow(t/300,-0.52);
}
double k112(double t) {
	return 7e-8*pow(t/300,-0.5);
}
double k113(double t) {
	return 1.6e-7*pow(t/300,-0.6);
}
double k114(double t) {
	return 4.03e-7*pow(t/300,-0.6);
}
double k115(double t) {
	return 7.68e-8*pow(t/300,-0.6);
}
double k116(double t) {
	return 7.75e-8*pow(t/300,-0.5);
}
double k117(double t) {
	return 1.95e-7*pow(t/300,-0.5);
}
double k118(double t) {
	return 2.0e-7*pow(t/300,-0.4);
}
double k119(double t) {
	return 6.3e-9*pow(t/300,-0.48);
}
double k120(double t) {
	return 3.05e-7*pow(t/300,-0.5);
}
double k121(double t) {
	return 3.9e-8*pow(t/300,-0.5);
}
double k122(double t) {
	return 8.6e-8*pow(t/300,-0.5);
}
double k123(double t) {
	return 1.08e-7*pow(t/300,-0.5);
}
double k124(double t) {
	return 6.02e-8*pow(t/300,-0.5);
}
double k125(double t) {
	return 2.58e-7*pow(t/300,-0.5);
}
double k126(double t) {
	return 5.6e-9*pow(t/300,-0.5);
}
double k127(double t) {
	return 1.95e-7*pow(t/300,-0.7);
}
double k128(double t) {
	return 2.75e-7*pow(t/300,-0.55);
}
double k129(double t) {
	return 2.76e-7*pow(t/300,-0.64);
}
double k130(double t) {
	return 2.4e-8*pow(t/300,-0.64);
}
double k131(double t) {
	return 1.1e-7*pow(t/300,-1.0);
}
double k132(double t) {
	return 1e-9;
}
double k133(double t) {
	return 1e-9;
}
double k134(double t) {
	return 1e-10;
}
double k135(double t) {
	return 5e-10;
}
double k136(double t) {
	return 1e-13;
}
double k137(double t) {
	return 5e-10;
}
double k138(double t) {
	return 5e-10;
}
double k139(double t) {
	return 7e-10;
}
double k140(double t) {
	return 5e-10;
}
double k141(double t) {
	return 1e-16;
}
double k142(double t) {
	return 2.25e-15;
}
double k143(double t) {
	return 1e-17;
}
double k144(double t) {
	return 1e-17;
}
double k145(double t) {
	return 4.36e-18*pow(t/300,0.35)*exp(-161.3/t);
}
double k146(double t) {
	if (t <= 300) return 2.1e-19;
	else return 3.09e-17*pow(t/300,0.33)*exp(-1629/t);
}
double k147(double t) {
	return 4.46e-16*pow(t,-0.5)*exp(-4.93/pow(t,0.6667));
}
double k148(double t) {
	return 4e-16*pow(t/300,-0.2);
}
double k149(double t) {
	if (t <= 300) return  2.5e-18;
	else return 3.14e-18*pow(t/300,-0.15)*exp(68/t);
}
double k150(double t) {
	return 1.5e-15;
}
double k151(double t) {
	return 9.9e-19*pow(t/300,-0.38);
}
double k152(double t) {
	return 4.9e-20*pow(t/300,1.58);
}
double k153(double t) {
	return 5.26e-18*pow(t/300,-5.22)*exp(-90/t);
}
double k154(double t) {
	if (t <= 300) return 1.32e-32*pow(t/300,-0.38);
	else return 1.32e-32*pow(t/300,-1.0);
}
double k155(double t) {
	return 2.8e-31*pow(t,-0.6);
}
double k156(double t) {
	return 6.9e-32*pow(t,-0.4);
}
double k157(double t) {
	if (t <= 5000) return 5.99e-33*pow(t/5000,-1.6);
	else return 5.99e-33*pow(t/5000,-0.64)*exp(5255/t);
}
double k158(double t) {
	if (t < 2000) return 6.16e-29*pow(t/300,-3.08);
	else return 2.14e-29*pow(t/300,-3.08)*exp(2114/t);
}
double k159(double t) {
	return 
}
double k160(double t) {
	return 
}
double k161(double t) {
	return 4.33e-32*pow(t/300,-1.0);
}
double k162(double t) {
	return 2.56e-31*pow(t/300,-2.0);
}
double k163(double t) {
	return 9.2e-34*pow(t/300,-1.0);
}
double k164(double t) {
	return 2.0e-11*pow(t/300,0.44);
}
double k165(double t, double temp_grain) {
	double fA = 1/(1.0 + 1e4*exp(-600/temp_grain));
	return 1/(3.0e-18*sqrt(t)*fA*(1.0 + 0.04*sqrt(t + temp_grain) + 0.002*t + 8e-6*t*t));
}
