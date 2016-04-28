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

double getxH2(double vec_pHx[]);
double getxH(double vec_pHx[]);

void derivs(double t, int nvar, double vec_nHx[], double vec_dnHxdt[]);
void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar);

double getTemp();
void setTemp(int t);

double getGrainTemp();
void setGrainTemp(int t);

double k1(double t);
double k2(double t);
double k3(double t);
double k4(double t);
double k5(double t);
double k6(double t);
double k7(double t);
double k8(double t);
double k9(double t, double nH);
double k10(double t, double nH2);
double k11(double t);
double k12(double t);
double k13(double t);
double k14(double t);
double k15(double t);
double k16(double t);
double k17(double t);
double k18(double t);
double k19(double t);
double k20(double t);
double k21(double t);
double k22(double t);
double k23(double t);
double k24(double t);
double k25(double t);
double k26(double t);
double k27(double t);
double k28(double t);
double k29(double t);
double k30(double t, double nH2);
double k31(double t);
double k32(double t);
double k33(double t);
double k34(double t);
double k35(double t);
double k36(double t);
double k37(double t);
double k38(double t);
double k39(double t);
double k40(double t);
double k41(double t);
double k42(double t);
double k43(double t);
double k44(double t);
double k45(double t);
double k46(double t);
double k47(double t);
double k48(double t);
double k49(double t);
double k50(double t);
double k51(double t);
double k52(double t);
double k53(double t);
double k54(double t);
double k55(double t);
double k56(double t);
double k57(double t);
double k58(double t);
double k59(double t);
double k60(double t);
double k61(double t);
double k62(double t);
double k63(double t);
double k64(double t);
double k65(double t);
double k66(double t);
double k67(double t);
double k68(double t);
double k69(double t);
double k70(double t);
double k71(double t);
double k72(double t);
double k73(double t);
double k74(double t);
double k75(double t);
double k76(double t);
double k77(double t);
double k78(double t);
double k79(double t);
double k80(double t);
double k81(double t);
double k82(double t);
double k83(double t);
double k84(double t);
double k85(double t);
double k86(double t);
double k87(double t);
double k88(double t);
double k89(double t);
double k90(double t);
double k91(double t);
double k92(double t);
double k93(double t);
double k94(double t);
double k95(double t);
double k96(double t);
double k97(double t);
double k98(double t);
double k99(double t);
double k100(double t);
double k101(double t);
double k102(double t);
double k103(double t);
double k104(double t);
double k105(double t);
double k106(double t);
double k107(double t);
double k108(double t);
double k109(double t);
double k110(double t);
double k111(double t);
double k112(double t);
double k113(double t);
double k114(double t);
double k115(double t);
double k116(double t);
double k117(double t);
double k118(double t);
double k119(double t);
double k120(double t);
double k121(double t);
double k122(double t);
double k123(double t);
double k124(double t);
double k125(double t);
double k126(double t);
double k127(double t);
double k128(double t);
double k129(double t);
double k130(double t);
double k131(double t);
double k132(double t);
double k133(double t);
double k134(double t);
double k135(double t);
double k136(double t);
double k137(double t);
double k138(double t);
double k139(double t);
double k140(double t);
double k141(double t);
double k142(double t);
double k143(double t);
double k144(double t);
double k145(double t);
double k146(double t);
double k147(double t);
double k148(double t);
double k149(double t);
double k150(double t);
double k151(double t);
double k152(double t);
double k153(double t);
double k154(double t);
double k155(double t);
double k156(double t);
double k157(double t);
double k158(double t);
double k159(double t);
double k160(double t);
double k161(double t);
double k162(double t);
double k163(double t);
double k164(double t);
double k165(double t, double temp_grain);


#endif