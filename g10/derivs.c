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
double k9(double t) {
	return 
}
double k10(double t) {
	return 
}
double k11(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-32.71396786 + 13.5365560*lt - 5.73932875*lt*lt + 1.56315498*lt*lt*lt - 2.87705600e−1*lt*lt*lt*lt + 3.48255977e−2*lt*lt*lt*lt*lt - 2.63197617e−3*lt*lt*lt*lt*lt*ln + 1.11954395e−4*lt*lt*lt*lt*lt*ln*ln - 2.03914985e−6*lt*lt*lt*lt*lt*ln*ln*ln);
}
double k12(double t) {
	return 
}
double k13(double t) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-1.801849334e1 + 2.36085220*lt - 2.82744300e−1*lt*lt + 1.62331664e−2*lt*lt*lt - 3.36501203e−2*lt*lt*lt*lt + 1.17832978e−2*lt*lt*lt*lt*lt - 1.65619470e−3*lt*lt*lt*lt*lt*lt + 1.06827520e−4*lt*lt*lt*lt*lt*lt*lt - 2.63128581e−6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k14(double t) {
	return 
}
double k15(double t) {
	return 
}
double k16(double t) {
	return 
}
double k17(double t) {
	return 
}
double k18(double t) {
	return 
}
double k19(double t) {
	return 
}
double k20(double t) {
	return 
}
double k21(double t) {
	return 
}
double k22(double t) {
	return 
}
double k23(double t) {
	return 
}
double k24(double t) {
	return 
}
double k25(double t) {
	return 
}
double k26(double t) {
	return 
}
double k27(double t) {
	return 
}
double k28(double t) {
	return 
}
double k29(double t) {
	return 
}
double k30(double t) {
	return 
}
double k31(double t) {
	return 
}
double k32(double t) {
	return 
}
double k33(double t) {
	return 
}
double k34(double t) {
	return 
}
double k35(double t) {
	return 
}
double k36(double t) {
	return 
}
double k37(double t) {
	return 
}
double k38(double t) {
	return 
}
double k39(double t) {
	return 
}
double k40(double t) {
	return 
}
double k41(double t) {
	return 
}
double k42(double t) {
	return 
}
double k43(double t) {
	return 
}
double k44(double t) {
	return 
}
double k45(double t) {
	return 
}
double k46(double t) {
	return 
}
double k47(double t) {
	return 
}
double k48(double t) {
	return 
}
double k49(double t) {
	return 
}
double k50(double t) {
	return 
}
double k51(double t) {
	return 
}
double k52(double t) {
	return 
}
double k53(double t) {
	return 
}
double k54(double t) {
	return 
}
double k55(double t) {
	return 
}
double k56(double t) {
	return 
}
double k57(double t) {
	return 
}
double k58(double t) {
	return 
}
double k59(double t) {
	return 
}
double k60(double t) {
	return 
}
double k61(double t) {
	return 
}
double k62(double t) {
	return 
}
double k63(double t) {
	return 
}
double k64(double t) {
	return 
}
double k65(double t) {
	return 
}
double k66(double t) {
	return 
}
double k67(double t) {
	return 
}
double k68(double t) {
	return 
}
double k69(double t) {
	return 
}
double k70(double t) {
	return 
}
double k71(double t) {
	return 
}
double k72(double t) {
	return 
}
double k73(double t) {
	return 
}
double k74(double t) {
	return 
}
double k75(double t) {
	return 
}
double k76(double t) {
	return 
}
double k77(double t) {
	return 
}
double k78(double t) {
	return 
}
double k79(double t) {
	return 
}
double k80(double t) {
	return 
}
double k81(double t) {
	return 
}
double k82(double t) {
	return 
}
double k83(double t) {
	return 
}
double k84(double t) {
	return 
}
double k85(double t) {
	return 
}
double k86(double t) {
	return 
}
double k87(double t) {
	return 
}
double k88(double t) {
	return 
}
double k89(double t) {
	return 
}
double k90(double t) {
	return 
}
double k91(double t) {
	return 
}
double k92(double t) {
	return 
}
double k93(double t) {
	return 
}
double k94(double t) {
	return 
}
double k95(double t) {
	return 
}
double k96(double t) {
	return 
}
double k97(double t) {
	return 
}
double k98(double t) {
	return 
}
double k99(double t) {
	return 
}
double k100(double t) {
	return 
}
double k101(double t) {
	return 
}
double k102(double t) {
	return 
}
double k103(double t) {
	return 
}
double k104(double t) {
	return 
}
double k105(double t) {
	return 
}
double k106(double t) {
	return 
}
double k107(double t) {
	return 
}
double k108(double t) {
	return 
}
double k109(double t) {
	return 
}
double k110(double t) {
	return 
}
double k111(double t) {
	return 
}
double k112(double t) {
	return 
}
double k113(double t) {
	return 
}
double k114(double t) {
	return 
}
double k115(double t) {
	return 
}
double k116(double t) {
	return 
}
double k117(double t) {
	return 
}
double k118(double t) {
	return 
}
double k119(double t) {
	return 
}
double k120(double t) {
	return 
}
double k121(double t) {
	return 
}
double k122(double t) {
	return 
}
double k123(double t) {
	return 
}
double k124(double t) {
	return 
}
double k125(double t) {
	return 
}
double k126(double t) {
	return 
}
double k127(double t) {
	return 
}
double k128(double t) {
	return 
}
double k129(double t) {
	return 
}
double k130(double t) {
	return 
}
double k131(double t) {
	return 
}
double k132(double t) {
	return 
}
double k133(double t) {
	return 
}
double k134(double t) {
	return 
}
double k135(double t) {
	return 
}
double k136(double t) {
	return 
}
double k137(double t) {
	return 
}
double k138(double t) {
	return 
}
double k139(double t) {
	return 
}
double k140(double t) {
	return 
}
double k141(double t) {
	return 
}
double k142(double t) {
	return 
}
double k143(double t) {
	return 
}
double k144(double t) {
	return 
}
double k145(double t) {
	return 
}
double k146(double t) {
	return 
}
double k147(double t) {
	return 
}
double k148(double t) {
	return 
}
double k149(double t) {
	return 
}
double k150(double t) {
	return 
}
double k151(double t) {
	return 
}
double k152(double t) {
	return 
}
double k153(double t) {
	return 
}
double k154(double t) {
	return 
}
double k155(double t) {
	return 
}
double k156(double t) {
	return 
}
double k157(double t) {
	return 
}
double k158(double t) {
	return 
}
double k159(double t) {
	return 
}
double k160(double t) {
	return 
}
double k161(double t) {
	return 
}
double k162(double t) {
	return 
}
double k163(double t) {
	return 
}
double k164(double t) {
	return 
}
double k165(double t) {
	return 
}
