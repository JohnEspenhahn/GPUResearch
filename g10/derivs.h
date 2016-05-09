#ifndef _DERIVS_H_
#define _DERIVS_H_

#include <math.h>
#include <float.h>
#include "../constants.h"
// #include <cmath.h>
#include "../jacobian.h"

#define ZONES 1
#define NDERV 27
#define NVAR (NDERV*ZONES)

#define N_tot 300.0  // cm^-3

#define xSi 1.5e-5
#define nSi (N_tot*xSi)

/*
#define xH_m 1e-6
#define xH2_p 1e-6
#define xH3_p 1e-6
#define xCH_p 1e-6
#define xCH2_p 1e-6
#define xOH_p 1e-6
#define xH2O_p 1e-6
#define xH3O_p 1e-6
#define xCO_p 1e-6
#define xHOC_p 1e-6
#define xO_m 1e-6
#define xC_m 1e-6
#define xO2_p 1e-6

#define nH_m (N_tot*xH_m)
#define nH2_p (N_tot*xH2_p)
#define nH3_p (N_tot*xH3_p)
#define nCH_p (N_tot*xCH_p)
#define nCH2_p (N_tot*xCH2_p)
#define nOH_p (N_tot*xOH_p)
#define nH2O_p (N_tot*xH2O_p)
#define nH3O_p (N_tot*xH3O_p)
#define nCO_p (N_tot*xCO_p)
#define nHOC_p (N_tot*xHOC_p)
#define nO_m (N_tot*xO_m)
#define nC_m (N_tot*xC_m)
#define nO2_p (N_tot*xO2_p)
*/

// Define total abundance for dissociated species
#define xSi_tot xSi
#define xC_tot  1.41e-4
#define xO_tot  3.16e-4
#define xHe_tot 0.1

#define xH_tot  (1.0 - xHe_tot - xC_tot - xO_tot - xSi_tot)

// Define total number density of dissociated species
#define N_H_tot (N_tot*xH_tot)
#define N_He_tot (N_tot*xHe_tot)
#define N_C_tot (N_tot*xC_tot)
#define N_O_tot (N_tot*xO_tot)

#define TEMP_INIT 50.0 // K
#define GRAIN_TEMP_INIT 50.0 // K

#define CDENSITY_H_tot 8e21 // cm^-2
#define SELF_SHIELDING (exp(-2e-21 * CDENSITY_H_tot) * 3.3e-11*1.7) // s^-1

#define COSMIC_RAY_RATE 1e-17*STEP_TIME // s^-1 try  try 1.8e-17

#define STEP_TIME 1 // s

double getnH(double[]);
double getnHe(double[]);
double getnC(double[]);
double getnO(double[]);
double getne(double[]);

double getxH2(double vec_pHx[]);
double getxH(double vec_pHx[]);

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getTemp();
void setTemp(int t);

double getGrainTemp();
void setGrainTemp(int t);

double k1(double, double[]);
double k2(double, double[]);
double k3(double, double[]);
double k4(double, double[]);
double k5(double, double[]);
double k6(double, double[]);
double k7(double, double[]);
double k8(double, double[]);
double k9(double, double[]);
double k10(double, double[]);
double k11(double, double[]);
double k12(double, double[]);
double k13(double, double[]);
double k14(double, double[]);
double k15(double, double[]);
double k16(double, double[]);
double k17(double, double[]);
double k18(double, double[]);
double k19(double, double[]);
double k20(double, double[]);
double k21(double, double[]);
double k22(double, double[]);
double k23(double, double[]);
double k24(double, double[]);
double k25(double, double[]);
double k26(double, double[]);
double k27(double, double[]);
double k28(double, double[]);
double k29(double, double[]);
double k30(double, double[]);
double k31(double, double[]);
double k32(double, double[]);
double k33(double, double[]);
double k34(double, double[]);
double k35(double, double[]);
double k36(double, double[]);
double k37(double, double[]);
double k38(double, double[]);
double k39(double, double[]);
double k40(double, double[]);
double k41(double, double[]);
double k42(double, double[]);
double k43(double, double[]);
double k44(double, double[]);
double k45(double, double[]);
double k46(double, double[]);
double k47(double, double[]);
double k48(double, double[]);
double k49(double, double[]);
double k50(double, double[]);
double k51(double, double[]);
double k52(double, double[]);
double k53(double, double[]);
double k54(double, double[]);
double k55(double, double[]);
double k56(double, double[]);
double k57(double, double[]);
double k58(double, double[]);
double k59(double, double[]);
double k60(double, double[]);
double k61(double, double[]);
double k62(double, double[]);
double k63(double, double[]);
double k64(double, double[]);
double k65(double, double[]);
double k66(double, double[]);
double k67(double, double[]);
double k68(double, double[]);
double k69(double, double[]);
double k70(double, double[]);
double k71(double, double[]);
double k72(double, double[]);
double k73(double, double[]);
double k74(double, double[]);
double k75(double, double[]);
double k76(double, double[]);
double k77(double, double[]);
double k78(double, double[]);
double k79(double, double[]);
double k80(double, double[]);
double k81(double, double[]);
double k82(double, double[]);
double k83(double, double[]);
double k84(double, double[]);
double k85(double, double[]);
double k86(double, double[]);
double k87(double, double[]);
double k88(double, double[]);
double k89(double, double[]);
double k90(double, double[]);
double k91(double, double[]);
double k92(double, double[]);
double k93(double, double[]);
double k94(double, double[]);
double k95(double, double[]);
double k96(double, double[]);
double k97(double, double[]);
double k98(double, double[]);
double k99(double, double[]);
double k100(double, double[]);
double k101(double, double[]);
double k102(double, double[]);
double k103(double, double[]);
double k104(double, double[]);
double k105(double, double[]);
double k106(double, double[]);
double k107(double, double[]);
double k108(double, double[]);
double k109(double, double[]);
double k110(double, double[]);
double k111(double, double[]);
double k112(double, double[]);
double k113(double, double[]);
double k114(double, double[]);
double k115(double, double[]);
double k116(double, double[]);
double k117(double, double[]);
double k118(double, double[]);
double k119(double, double[]);
double k120(double, double[]);
double k121(double, double[]);
double k122(double, double[]);
double k123(double, double[]);
double k124(double, double[]);
double k125(double, double[]);
double k126(double, double[]);
double k127(double, double[]);
double k128(double, double[]);
double k129(double, double[]);
double k130(double, double[]);
double k131(double, double[]);
double k132(double, double[]);
double k133(double, double[]);
double k134(double, double[]);
double k135(double, double[]);
double k136(double, double[]);
double k137(double, double[]);
double k138(double, double[]);
double k139(double, double[]);
double k140(double, double[]);
double k141(double, double[]);
double k142(double, double[]);
double k143(double, double[]);
double k144(double, double[]);
double k145(double, double[]);
double k146(double, double[]);
double k147(double, double[]);
double k148(double, double[]);
double k149(double, double[]);
double k150(double, double[]);
double k151(double, double[]);
double k152(double, double[]);
double k153(double, double[]);
double k154(double, double[]);
double k155(double, double[]);
double k156(double, double[]);
double k157(double, double[]);
double k158(double, double[]);
double k159(double, double[]);
double k160(double, double[]);
double k161(double, double[]);
double k162(double, double[]);
double k163(double, double[]);
double k164(double, double[]);
double k165(double t, double temp_grain);


#endif