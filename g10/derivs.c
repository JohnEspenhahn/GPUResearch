#include "derivs.h"

int T = TEMP_INIT, GRAIN_TEMP = GRAIN_TEMP_INIT;

double getnH(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27];
	
	return N_H_tot - nH_p - 2*nH2 - nOH - 2*nH2O - nHCO_p - nCH - 2*nCH2 - 3*nCH3_p - nH_m - 2*nH2_p - 3*nH3_p - nCH_p - 2*nCH2_p - nOH_p - 2*nH2O_p - 3*nH3O_p - nHOC_p;
}

double getnHe(double sub_vec_nHx[]) {
	double nHe_p = sub_vec_nHx[3];
	return N_He_tot - nHe_p;
}

double getnC(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27];
	
	return N_C_tot - nC_p - nCO - 2*nC2 - nHCO_p - nCH - nCH2 - nCH3_p - nCH_p - nCH2_p - nCO_p - nHOC_p - nC_m;
}

double getnO(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27];
	
	return N_O_tot - nO_p - nOH - nH2O - nCO - 2*nO2 - nHCO_p - nOH_p - nH2O_p - nH3O_p - nCO_p - nHOC_p - nO_m - 2*nO2_p;
}

double getne(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27];
	
	return -nH_m + nH2_p + nH3_p + nCH_p + nCH2_p + nOH_p + nH2O_p + nH3O_p + nCO_p + nHOC_p - nO_m - nC_m + nO2_p 
				+ nH_p + nHe_p + nC_p + nO_p + nHCO_p + nCH3_p;
}

double getxH(double sub_vec_nHx[]) {
	return getnH(sub_vec_nHx) / N_tot;
}

double getxH2(double sub_vec_nHx[]) {
	double nH2 = sub_vec_nHx[2];
	return (2*nH2) / N_tot;
}

double getpH_p(double nH_p) { return nH_p * (mu_h * M_h); }

double getpH2(double nH2) { return nH2 * (mu_h2 * M_h); }

double getpHe_p(double nHe_p) { return nHe_p * (mu_he * M_h); }

double getpC_p(double nC_p) { return nC_p * (mu_c * M_h); }

double getpO_p(double nO_p) { return nO_p * (mu_o * M_h); }

double getpOH(double nOH) { return nOH * (mu_oh * M_h); }

double getpH2O(double nH2O) { return nH2O * (mu_h2o * M_h); }

double getpCO(double nCO) { return nCO * (mu_co * M_h); }

double getpC2(double nC2) { return nC2 * (mu_c2 * M_h); }

double getpO2(double nO2) { return nO2 * (mu_o2 * M_h); }

double getpHCO_p(double nHCO_p) { return nHCO_p * (mu_hco * M_h); }

double getpCH(double nCH) { return nCH * (mu_ch * M_h); }

double getpCH2(double nCH2) { return nCH2 * (mu_ch2 * M_h); }

double getpCH3_p(double nCH3_p) { return nCH3_p * (mu_ch3 * M_h); }

double dnO_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k23(t)*ne*nO + k25(t)*nH_p*nO + k26(t)*nHe_p*nO + k81(t)*nO2*nC_p + k95(t)*nOH*nHe_p + k102(t)*nO2*nHe_p + k105(t)*nCO*nHe_p - k21(t)*nO_p*ne - k24(t)*nO_p*nH - k68(t)*nO_p*nC2 - k69(t)*nO_p*nH2 - k160(t)*nO_p*nC*nM;
}

double dnCO_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k61(t)*nCH_p*nO + k68(t)*nO_p*nC2 + k73(t)*nC_p*nOH + k80(t)*nO2*nC_p + k83(t)*nO2_p*nC + k149(t)*nC_p*nO + k159(t)*nC_p*nM*nO + k160(t)*nO_p*nC*nM - k106(t)*nCO_p*nH - k128(t)*nCO_p*ne;
}

double dnO2(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k47(t)*nOH*nO + k103(t)*nO2_p*nC + k152(t)*nO*nO + k163(t)*nM*nO*nO - k50(t)*nO2*nH - k51(t)*nO2*nH2 - k52(t)*nO2*nC - k80(t)*nO2*nC_p - k81(t)*nO2*nC_p - k82(t)*nO2*nCH2_p - k100(t)*nO2*nH_p - k101(t)*nO2*nHe_p - k102(t)*nO2*nHe_p;
}

double dnCH_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k56(t)*nC*nH2_p + k57(t)*nC*nH3_p + k58(t)*nC_p*nH2 + k62(t)*nH_p*nCH2 + k63(t)*nH*nCH2_p + k86(t)*nHCO_p*nC + k90(t)*nCH*nH_p + k147(t)*nC_p*nH - k59(t)*nH*nCH_p - k60(t)*nCH_p*nH2 - k61(t)*nCH_p*nO - k112(t)*ne*nCH_p;
}

double dnCH3_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k64(t)*nH2*nCH2_p - k66(t)*nH*nCH3_p - k67(t)*nCH3_p*nO - k116(t)*ne*nCH3_p - k117(t)*ne*nCH3_p - k118(t)*ne*nCH3_p;
}

double dnC2(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k37(t)*nC*nCH + k145(t)*nC*nC + k157(t)*nC*nC*nM - k42(t)*nC2*nO - k68(t)*nO_p*nC2 - k93(t)*nHe_p*nC2;
}

double dnCH2_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k60(t)*nCH_p*nH2 + k66(t)*nH*nCH3_p + k91(t)*nH_p*nCH2 + k148(t)*nC_p*nH2 - k63(t)*nH*nCH2_p - k64(t)*nH2*nCH2_p - k65(t)*nCH2_p*nO - k82(t)*nO2*nCH2_p - k113(t)*ne*nCH2_p - k114(t)*ne*nCH2_p - k115(t)*ne*nCH2_p;
}

double dnO2_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k100(t)*nO2*nH_p + k101(t)*nO2*nHe_p - k83(t)*nO2_p*nC - k103(t)*nO2_p*nC - k127(t)*nO2_p*ne;
}

double dnHOC_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k78(t)*nH2O*nC_p + k84(t)*nCO*nH3_p - k32(t)*nHOC_p*nH2 - k33(t)*nHOC_p*nCO - k131(t)*ne*nHOC_p;
}

double dnH2O_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k72(t)*nOH*nH3_p + k74(t)*nH2*nOH_p + k96(t)*nH2O*nH_p + k99(t)*nH2O*nHe_p - k75(t)*nH2O_p*nH2 - k120(t)*ne*nH2O_p - k121(t)*ne*nH2O_p - k122(t)*ne*nH2O_p;
}

double dnOH(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k43(t)*nH2*nO + k49(t)*nH2O*nH + k50(t)*nO2*nH + k51(t)*nO2*nH2 + k53(t)*nH*nCO + k82(t)*nO2*nCH2_p + k97(t)*nH2O*nHe_p + k122(t)*ne*nH2O_p + k124(t)*ne*nH3O_p + k125(t)*ne*nH3O_p + k130(t)*nHCO_p*ne + k133(t)*nH_m*nO + k138(t)*nH*nO_m + k151(t)*nH*nO + k161(t)*nH*nM*nO - k31(t)*nH*nOH - k44(t)*nH*nOH - k45(t)*nOH*nH2 - k46(t)*nC*nOH - k47(t)*nOH*nO - k48(t)*nOH*nOH - k72(t)*nOH*nH3_p - k73(t)*nC_p*nOH - k94(t)*nH_p*nOH - k95(t)*nOH*nHe_p - k134(t)*nH_m*nOH - k153(t)*nH*nOH - k162(t)*nH*nOH*nM;
}

double dnH2O(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k45(t)*nOH*nH2 + k48(t)*nOH*nOH + k123(t)*ne*nH3O_p + k134(t)*nH_m*nOH + k139(t)*nH2*nO_m + k153(t)*nH*nOH + k162(t)*nH*nOH*nM - k49(t)*nH2O*nH - k76(t)*nH2O*nH3_p - k77(t)*nH2O*nC_p - k78(t)*nH2O*nC_p - k87(t)*nHCO_p*nH2O - k96(t)*nH2O*nH_p - k97(t)*nH2O*nHe_p - k98(t)*nH2O*nHe_p - k99(t)*nH2O*nHe_p;
}

double dnCH(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k34(t)*nC*nH2 + k39(t)*nH*nCH2 + k113(t)*ne*nCH2_p + k117(t)*ne*nCH3_p + k118(t)*ne*nCH3_p + k132(t)*nC*nH_m + k135(t)*nC_m*nH + k143(t)*nC*nH - k35(t)*nCH*nH - k36(t)*nCH*nH2 - k37(t)*nC*nCH - k38(t)*nCH*nO - k90(t)*nCH*nH_p - k164(t)*nCH*nO;
}

double dnH2(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T
			 , grain_temp = GRAIN_TEMP;

	return k2(t)*nH_m*nH + k4(t)*nH*nH2_p + k35(t)*nCH*nH + k39(t)*nH*nCH2 + k41(t)*nCH2*nO + k44(t)*nH*nOH + k49(t)*nH2O*nH + k55(t)*nH*nH3_p + k57(t)*nC*nH3_p + k59(t)*nH*nCH_p + k62(t)*nH_p*nCH2 + k63(t)*nH*nCH2_p + k66(t)*nH*nCH3_p + k67(t)*nCH3_p*nO + k71(t)*nH3_p*nO + k72(t)*nOH*nH3_p + k76(t)*nH2O*nH3_p + k79(t)*nC*nH3O_p + k84(t)*nCO*nH3_p + k85(t)*nCO*nH3_p + k92(t)*nCH2*nHe_p + k110(t)*ne*nH3_p + k115(t)*ne*nCH2_p + k117(t)*ne*nCH3_p + k121(t)*ne*nH2O_p + k124(t)*ne*nH3O_p + k126(t)*ne*nH3O_p + k154(t)*nH*nH*nH + k155(t)*nH*nH*nH2 + k156(t)*nH*nH*nHe + k165(t,grain_temp)*nH*nH - k7(t)*nH_p*nH2 - k8(t)*ne*nH2 - k9(t,nH)*nH*nH2 - k10(t,nH2)*nH2*nH2 - k30(t,nH2)*nH2*nHe - k34(t)*nC*nH2 - k36(t)*nCH*nH2 - k43(t)*nH2*nO - k45(t)*nOH*nH2 - k51(t)*nO2*nH2 - k54(t)*nH2_p*nH2 - k58(t)*nC_p*nH2 - k60(t)*nCH_p*nH2 - k64(t)*nH2*nCH2_p - k69(t)*nO_p*nH2 - k74(t)*nH2*nOH_p - k75(t)*nH2O_p*nH2 - k88(t)*nH2*nHe_p - k89(t)*nH2*nHe_p - k136(t)*nC_m*nH2 - k139(t)*nH2*nO_m - k141(t)*nH_p*nH2 - k144(t)*nC*nH2 - k148(t)*nC_p*nH2;
}

double dnHe_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k16(t)*ne*nHe + k19(t)*nH_p*nHe - k17(t)*ne*nHe_p - k18(t)*nH*nHe_p - k26(t)*nHe_p*nO - k29(t)*nC*nHe_p - k88(t)*nH2*nHe_p - k89(t)*nH2*nHe_p - k92(t)*nCH2*nHe_p - k93(t)*nHe_p*nC2 - k95(t)*nOH*nHe_p - k97(t)*nH2O*nHe_p - k98(t)*nH2O*nHe_p - k99(t)*nH2O*nHe_p - k101(t)*nO2*nHe_p - k102(t)*nO2*nHe_p - k104(t)*nCO*nHe_p - k105(t)*nCO*nHe_p - k109(t)*nH_m*nHe_p;
}

double dnCO(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k38(t)*nCH*nO + k40(t)*nCH2*nO + k41(t)*nCH2*nO + k42(t)*nC2*nO + k46(t)*nC*nOH + k52(t)*nO2*nC + k81(t)*nO2*nC_p + k86(t)*nHCO_p*nC + k87(t)*nHCO_p*nH2O + k106(t)*nCO_p*nH + k129(t)*nHCO_p*ne + k131(t)*ne*nHOC_p + k137(t)*nC_m*nO + k140(t)*nC*nO_m + k146(t)*nC*nO + k158(t)*nC*nM*nO - k53(t)*nH*nCO - k84(t)*nCO*nH3_p - k85(t)*nCO*nH3_p - k104(t)*nCO*nHe_p - k105(t)*nCO*nHe_p;
}

double dnCH2(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k36(t)*nCH*nH2 + k116(t)*ne*nCH3_p + k136(t)*nC_m*nH2 + k144(t)*nC*nH2 - k39(t)*nH*nCH2 - k40(t)*nCH2*nO - k41(t)*nCH2*nO - k62(t)*nH_p*nCH2 - k91(t)*nH_p*nCH2 - k92(t)*nCH2*nHe_p;
}

double dnOH_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k69(t)*nO_p*nH2 + k70(t)*nH2_p*nO + k71(t)*nH3_p*nO + k94(t)*nH_p*nOH + k98(t)*nH2O*nHe_p - k74(t)*nH2*nOH_p - k119(t)*ne*nOH_p;
}

double dnC_m(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k142(t)*nC*ne - k107(t)*nC_m*nH_p - k135(t)*nC_m*nH - k136(t)*nC_m*nH2 - k137(t)*nC_m*nO;
}

double dnHCO_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k32(t)*nHOC_p*nH2 + k33(t)*nHOC_p*nCO + k65(t)*nCH2_p*nO + k67(t)*nCH3_p*nO + k77(t)*nH2O*nC_p + k79(t)*nC*nH3O_p + k82(t)*nO2*nCH2_p + k85(t)*nCO*nH3_p + k164(t)*nCH*nO - k86(t)*nHCO_p*nC - k87(t)*nHCO_p*nH2O - k129(t)*nHCO_p*ne - k130(t)*nHCO_p*ne;
}

double dnC_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k22(t)*nC*ne + k27(t)*nC*nH_p + k29(t)*nC*nHe_p + k59(t)*nH*nCH_p + k92(t)*nCH2*nHe_p + k93(t)*nHe_p*nC2 + k103(t)*nO2_p*nC + k104(t)*nCO*nHe_p - k20(t)*nC_p*ne - k28(t)*nC_p*nH - k58(t)*nC_p*nH2 - k73(t)*nC_p*nOH - k77(t)*nH2O*nC_p - k78(t)*nH2O*nC_p - k80(t)*nO2*nC_p - k81(t)*nO2*nC_p - k147(t)*nC_p*nH - k148(t)*nC_p*nH2 - k149(t)*nC_p*nO - k159(t)*nC_p*nM*nO;
}

double dnH_m(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k1(t)*ne*nH - k2(t)*nH_m*nH - k5(t)*nH_m*nH_p - k13(t)*ne*nH_m - k14(t)*nH_m*nH - k15(t)*nH_m*nH_p - k109(t)*nH_m*nHe_p - k132(t)*nC*nH_m - k133(t)*nH_m*nO - k134(t)*nH_m*nOH;
}

double dnH_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k4(t)*nH*nH2_p + k11(t)*ne*nH + k18(t)*nH*nHe_p + k24(t)*nO_p*nH + k28(t)*nC_p*nH + k89(t)*nH2*nHe_p + k97(t)*nH2O*nHe_p + k106(t)*nCO_p*nH - k3(t)*nH*nH_p - k5(t)*nH_m*nH_p - k7(t)*nH_p*nH2 - k12(t)*ne*nH_p - k15(t)*nH_m*nH_p - k19(t)*nH_p*nHe - k25(t)*nH_p*nO - k27(t)*nC*nH_p - k62(t)*nH_p*nCH2 - k90(t)*nCH*nH_p - k91(t)*nH_p*nCH2 - k94(t)*nH_p*nOH - k96(t)*nH2O*nH_p - k100(t)*nO2*nH_p - k107(t)*nC_m*nH_p - k108(t)*nH_p*nO_m - k141(t)*nH_p*nH2;
}

double dnH2_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k3(t)*nH*nH_p + k7(t)*nH_p*nH2 + k15(t)*nH_m*nH_p + k55(t)*nH*nH3_p + k88(t)*nH2*nHe_p - k4(t)*nH*nH2_p - k6(t)*ne*nH2_p - k54(t)*nH2_p*nH2 - k56(t)*nC*nH2_p - k70(t)*nH2_p*nO;
}

double dnH3O_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k75(t)*nH2O_p*nH2 + k76(t)*nH2O*nH3_p + k87(t)*nHCO_p*nH2O - k79(t)*nC*nH3O_p - k123(t)*ne*nH3O_p - k124(t)*ne*nH3O_p - k125(t)*ne*nH3O_p - k126(t)*ne*nH3O_p;
}

double dnO_m(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k150(t)*ne*nO - k108(t)*nH_p*nO_m - k138(t)*nH*nO_m - k139(t)*nH2*nO_m - k140(t)*nC*nO_m;
}

double dnH3_p(double sub_vec_nHx[]) {
	double nH_p = sub_vec_nHx[1], nH2  = sub_vec_nHx[2], nHe_p = sub_vec_nHx[3]
			 , nC_p = sub_vec_nHx[4], nO_p = sub_vec_nHx[5], nOH   = sub_vec_nHx[6]
			 , nH2O = sub_vec_nHx[7], nCO  = sub_vec_nHx[8], nC2   = sub_vec_nHx[9]
			 , nO2  = sub_vec_nHx[10], nHCO_p = sub_vec_nHx[11], nCH = sub_vec_nHx[12]
			 , nCH2 = sub_vec_nHx[13], nCH3_p = sub_vec_nHx[14], nH_m = sub_vec_nHx[15]
			 , nH2_p = sub_vec_nHx[16], nH3_p = sub_vec_nHx[17], nCH_p = sub_vec_nHx[18]
			 , nCH2_p = sub_vec_nHx[19], nOH_p = sub_vec_nHx[20], nH2O_p = sub_vec_nHx[21]
			 , nH3O_p = sub_vec_nHx[22], nCO_p = sub_vec_nHx[23], nHOC_p = sub_vec_nHx[24]
			 , nO_m = sub_vec_nHx[25], nC_m = sub_vec_nHx[26], nO2_p = sub_vec_nHx[27]
			 , nH = getnH(sub_vec_nHx)
			 , nC = getnC(sub_vec_nHx)
			 , ne = getne(sub_vec_nHx)
			 , nO = getnO(sub_vec_nHx)
			 , nHe = getnHe(sub_vec_nHx)
			 , nM = nC + nO + nSi
			 , t = T;

	return k54(t)*nH2_p*nH2 + k141(t)*nH_p*nH2 - k55(t)*nH*nH3_p - k57(t)*nC*nH3_p - k71(t)*nH3_p*nO - k72(t)*nOH*nH3_p - k76(t)*nH2O*nH3_p - k84(t)*nCO*nH3_p - k85(t)*nCO*nH3_p - k110(t)*ne*nH3_p - k111(t)*ne*nH3_p;
}

// run a step of the derivatives
void derivs(double x, int nvar, double vec_nHx[], double vec_dnHxdt[]) {
	for (int i = 1; i <= nvar; i += NDERV) {
		double *sub_vec_nHx = vec_nHx + (i-1);
		
		// pH_p
		vec_dnHxdt[i] = dnH_p(sub_vec_nHx);
		
		// pH2
		vec_dnHxdt[i+1] = dnH2(sub_vec_nHx);
		
		// pHe_p
		vec_dnHxdt[i+2] = dnHe_p(sub_vec_nHx);
		
		// pC_p
		vec_dnHxdt[i+3] = dnC_p(sub_vec_nHx);
		
		// pO_p
		vec_dnHxdt[i+4] = dnO_p(sub_vec_nHx);
		
		// pOH
		vec_dnHxdt[i+5] = dnOH(sub_vec_nHx);
		
		// pH2O
		vec_dnHxdt[i+6] = dnH2O(sub_vec_nHx);
		
		// pCO
		vec_dnHxdt[i+7] = dnCO(sub_vec_nHx);
		
		// pC2
		vec_dnHxdt[i+8] = dnC2(sub_vec_nHx);
		
		// pO2
		vec_dnHxdt[i+9] = dnO2(sub_vec_nHx);
		
		// pHCO_p
		vec_dnHxdt[i+10] = dnHCO_p(sub_vec_nHx);
		
		// pCH
		vec_dnHxdt[i+11] = dnCH(sub_vec_nHx);
		
		// pCH2
		vec_dnHxdt[i+12] = dnCH2(sub_vec_nHx);
		
		// pCH3_p
		vec_dnHxdt[i+13] = dnCH3_p(sub_vec_nHx);
		
		// H_m 
		vec_dnHxdt[i+14] = dnH_m(sub_vec_nHx);
		// H2_p 
		vec_dnHxdt[i+15] = dnH2_p(sub_vec_nHx);
		// H3_p 
		vec_dnHxdt[i+16] = dnH3_p(sub_vec_nHx);
		// CH_p 
		vec_dnHxdt[i+17] = dnCH_p(sub_vec_nHx);
		// CH2_p 
		vec_dnHxdt[i+18] = dnCH2_p(sub_vec_nHx);
		// OH_p 
		vec_dnHxdt[i+19] = dnOH_p(sub_vec_nHx);
		// H2O_p 
		vec_dnHxdt[i+20] = dnH2O_p(sub_vec_nHx);
		// H3O_p 
		vec_dnHxdt[i+21] = dnH3O_p(sub_vec_nHx);
		// CO_p 
		vec_dnHxdt[i+22] = dnCO_p(sub_vec_nHx);
		// HOC_p 
		vec_dnHxdt[i+23] = dnHOC_p(sub_vec_nHx);
		// O_m 
		vec_dnHxdt[i+24] = dnO_m(sub_vec_nHx);
		// C_m 
		vec_dnHxdt[i+25] = dnC_m(sub_vec_nHx);
		// O2_p 
		vec_dnHxdt[i+26] = dnO2_p(sub_vec_nHx);
	}
}

double (*derivs_arr[])(double[]) = { 
	&dnH_p, &dnH2, &dnHe_p, &dnC_p, &dnO_p, &dnOH, &dnH2O, &dnCO, &dnC2, &dnO2, &dnHCO_p, &dnCH, &dnCH2, &dnCH3_p
	, &dnH_m, &dnH2_p, &dnH3_p, &dnCH_p, &dnCH2_p, &dnOH_p, &dnH2O_p, &dnH3O_p, &dnCO_p, &dnHOC_p, &dnO_m, &dnC_m, &dnO2_p
};

void jacobn(double x, double vec_nHx[], double dfdx[], double **dfdy, int nvar)
{
	for (int i = 1; i <= nvar; i += NDERV) {
		// copy_vector(vec_nHx + (i-1), y_temp, 1, nvar);
		double *y_temp = vec_nHx + (i-1);		
		jacobian(derivs_arr, NDERV, y_temp, 1e-7, dfdy, nvar);
	}
}

double k1(double t, double sub_vec_nHx[]) {
	double lt = log10(t);
	if (t <= 6000) return pow(10, -17.845 + 0.762*lt + 0.1523*lt*lt - 0.03274*lt*lt*lt);
	else return pow(10, -16.420 + 0.1998*lt*lt - 5.447e-3*lt*lt*lt*lt + 4.0415e-5*lt*lt*lt*lt*lt*lt);
}
double k2(double t, double sub_vec_nHx[]) {
	if (t <= 300) return 1.5e-9;
	else return 4.0e-9*pow(t, -0.17);
}
double k3(double t, double sub_vec_nHx[]) {
	double lt = log10(t);
	return pow(10, -19.38 - 1.523*lt + 1.118*lt*lt - 0.1269*lt*lt*lt);
}
double k4(double t, double sub_vec_nHx[]) {
	return 6.4e-10;
}
double k5(double t, double sub_vec_nHx[]) {
	return 2.4e-6*pow(t, -0.5)*(1.0 + t/20000.0);
}
double k6(double t, double sub_vec_nHx[]) {
	if (t <= 617) return 1.0e-8;
	else return 1.32e-6*pow(t, -0.76);
}
double k7(double t, double sub_vec_nHx[]) {
	double lt = log(t); // ln
	double v = exp(-21237.15/t)*(-3.3232183e-7 + 3.3735382e-7*lt - 1.4491368e-7*lt*lt + 3.4172805e-8*lt*lt*lt - 4.7813720E-9*lt*lt*lt*lt + 3.9731542e-10*lt*lt*lt*lt*lt - 1.8171411e-11*lt*lt*lt*lt*lt*lt + 3.5311932e-13*lt*lt*lt*lt*lt*lt*lt);
	if (v < 0) return 0;
	else return v;
}
double k8(double t, double sub_vec_nHx[]) {
	return 3.73e-9*pow(t,0.1121)*exp(-99430.0/t);
}
double k9(double t, double sub_vec_nHx[]) {
	double kl = 6.67e-12*sqrt(t)*exp(-(1+63590/t)),
		   kh = 3.52e-9*exp(-43900/t),
		   logt = log10(t/10000),
		   nH = getnH(sub_vec_nHx),
		   ncr = nH / pow(10, 3.0 - 0.416*logt - 0.327*logt*logt);
		   
	if (kl < 1e-10) return 0;
	else return kh*pow(kh/kl,-1/(ncr+1));
}
double k10(double t, double sub_vec_nHx[]) {
	double kl = (5.996e-30*pow(t,4.1881)*exp(-54657.4/t))/pow(1.0 + 6.761e-6 * t, 5.6881),
		   kh = 1.3e-9*exp(-53300/t),
		   logt = log10(t/10000),
		   nH2 = sub_vec_nHx[2],
		   ncr = nH2 / pow(10, 4.845 - 1.3*logt + 1.62*logt*logt);
		   
	if (kl < 1e-10) return 0;
	return kh*pow(kh/kl,-1/(ncr+1));
}
double k11(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-32.71396786 + 13.5365560*lt - 5.73932875*lt*lt + 1.56315498*lt*lt*lt - 2.87705600e-1*lt*lt*lt*lt + 3.48255977e-2*lt*lt*lt*lt*lt - 2.63197617e-3*lt*lt*lt*lt*lt*lt + 1.11954395e-4*lt*lt*lt*lt*lt*lt*lt - 2.03914985e-6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k12(double t, double sub_vec_nHx[]) {
	return 1.269e-13*pow(315614/t, 1.503)*pow(1+pow(6046255/t,0.470),-1.923);
}
double k13(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-1.801849334e1 + 2.36085220*lt - 2.82744300e-1*lt*lt + 1.62331664e-2*lt*lt*lt - 3.36501203e-2*lt*lt*lt*lt + 1.17832978e-2*lt*lt*lt*lt*lt - 1.65619470e-3*lt*lt*lt*lt*lt*lt + 1.06827520e-4*lt*lt*lt*lt*lt*lt*lt - 2.63128581e-6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k14(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	if (te <= 0.1) return 2.5634e-9*pow(te,1.78186);
	else return exp(-2.0372609e1 + 1.13944933e0*lt - 1.4210135e-1*lt*lt + 8.4644554e-3*lt*lt*lt - 1.4327641e-3*lt*lt*lt*lt + 2.0122503e-4*lt*lt*lt*lt*lt + 8.6639632e-5*lt*lt*lt*lt*lt*lt - 2.5850097e-5*lt*lt*lt*lt*lt*lt*lt + 2.4555012e-6*lt*lt*lt*lt*lt*lt*lt*lt - 8.0683825e-8*lt*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k15(double t, double sub_vec_nHx[]) {
	if (t <= 8000) return 6.9e-9*pow(t,-0.35);
	else return 9.6e-7*pow(t,-0.90);
}
double k16(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double lt = log(te); // ln
	return exp(-4.409864886e1 + 2.391596563e1*lt - 1.07532302e1*lt*lt + 3.05803875e0*lt*lt*lt - 5.6851189e-1*lt*lt*lt*lt + 6.79539123e-2*lt*lt*lt*lt*lt - 5.0090561e-3*lt*lt*lt*lt*lt*lt + 2.06723616e-4*lt*lt*lt*lt*lt*lt*lt - 3.64916141e-6*lt*lt*lt*lt*lt*lt*lt*lt);
}
double k17(double t, double sub_vec_nHx[]) {
	double lt = log10(t); // log
	double di = 0; // TODO dielectric
	return 1e-11*pow(t,-0.5)*(11.19 - 1.676*lt - 0.2852*lt*lt + 0.04433*lt*lt*lt) + di;
}
double k18(double t, double sub_vec_nHx[]) {
	return 1.25e-15*pow(t/300,0.25);
}
double k19(double t, double sub_vec_nHx[]) {
	if (t <= 10000) return 1.26e-9*pow(t,-0.75)*exp(-127500/t);
	else return 4.0e-37*pow(t,4.74);
}
double k20(double t, double sub_vec_nHx[]) {
	if (t <= 7950) return 4.67e-12*pow(t/300,-0.6);
	else if (t <= 21140) return 1.23e-17*pow(t/300,2.49)*exp(21845.6/t);
	else return 9.62e-8*pow(t/300,-1.37)*exp(-115786.2/t);
}
double k21(double t, double sub_vec_nHx[]) {
	if (t <= 400) return 1.30e-10*pow(t,-0.64);
	else return 1.41e-10*pow(t,-0.66) + 7.4e-4*pow(t,-1.5)*exp(-175000/t)*(1.0 + 0.062*exp(-145000/t));
}
double k22(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double u = 11.26/te;
	return (6.85e-8/(0.193 + u)) * pow(u,0.25) * exp(-u);
}
double k23(double t, double sub_vec_nHx[]) {
	double te = t * 8.621738e-5;
	double u = 13.6/te;
	return (3.59e-8/(0.073 + u)) * pow(u,0.34) * exp(-u);
}
double k24(double t, double sub_vec_nHx[]) {
	return 4.99e-11 * pow(t,0.405) + 7.54e-10*pow(t,-0.458);
}
double k25(double t, double sub_vec_nHx[]) {
	return (1.08e-11*pow(t,0.517) + 4.00e-10*pow(t,0.00669)) * exp(-227/t);
}
double k26(double t, double sub_vec_nHx[]) {
	return 4.991e-15*pow(t/10000,0.3794)*exp(-t/1121000) + 2.78e-15*pow(t/10000,-0.2163)*exp(t/815800);
}
double k27(double t, double sub_vec_nHx[]) {
	return 3.9e-16*pow(t,0.213);
}
double k28(double t, double sub_vec_nHx[]) {
	return 6.08e-14*pow(t/10000,1.96)*exp(-170000/t);
}
double k29(double t, double sub_vec_nHx[]) {
	if (t <= 200) return 8.58e-17*pow(t,0.757);
	else if (t <= 2000) return 3.25e-17*pow(t,0.968);
	else return 2.77e-19*pow(t,1.597);
}
double k30(double t, double sub_vec_nHx[]) {
	double logt = log10(t),
		   kl = pow(10, -27.029 + 3.801*logt - 29487/t),
		   kh = pow(10, -2.729 - 1.75*logt - 23474/t),
		   nH2 = sub_vec_nHx[2],
		   ncr = pow(nH2 / pow(10, 5.0792*(1.0 - 1.23e-5*(t-2000))), 1.07);
		   
	if (kl < 1e-10) return 0;
	else return kh*pow(kh/kl,-1/(ncr+1));
}
double k31(double t, double sub_vec_nHx[]) {
	return 6.0e-9*exp(-50900/t);
}
double k32(double t, double sub_vec_nHx[]) {
	return 3.8e-10;
}
double k33(double t, double sub_vec_nHx[]) {
	return 4.0e-10;
}
double k34(double t, double sub_vec_nHx[]) {
	return 6.64e-10*exp(-11700/t);
}
double k35(double t, double sub_vec_nHx[]) {
	return 1.31e-10*exp(-80/t);
}
double k36(double t, double sub_vec_nHx[]) {
	return 5.46e-10*exp(-1943/t);
}
double k37(double t, double sub_vec_nHx[]) {
	return 6.59e-11;
}
double k38(double t, double sub_vec_nHx[]) {
	if (t <= 2000) return 6.6e-11;
	else return 1.02e-10*exp(-914/t);
}
double k39(double t, double sub_vec_nHx[]) {
	return 6.64e-11;
}
double k40(double t, double sub_vec_nHx[]) {
	return 1.33e-10;
}
double k41(double t, double sub_vec_nHx[]) {
	return 8.0e-11;
}
double k42(double t, double sub_vec_nHx[]) {
	if (t <= 300) return 5.0e-11*sqrt(t/500);
	else return 5.0e-11*pow(t/300,0.757);
}
double k43(double t, double sub_vec_nHx[]) {
	return 3.14e-13*pow(t/300,2.7)*exp(-3150/t);
}
double k44(double t, double sub_vec_nHx[]) {
	return 6.99e-14*pow(t/300,2.8)*exp(-1950/t);
}
double k45(double t, double sub_vec_nHx[]) {
	return 2.05e-12*pow(t/300,1.52)*exp(-1736/t);
}
double k46(double t, double sub_vec_nHx[]) {
	return 1e-10;
}
double k47(double t, double sub_vec_nHx[]) {
	if (t <= 261) return 3.5e-11;
	else return 1.77e-11*exp(178/t);
}
double k48(double t, double sub_vec_nHx[]) {
	return 1.65e-12*pow(t/300,1.14)*exp(-50/t);
}
double k49(double t, double sub_vec_nHx[]) {
	return 1.59e-11*pow(t/300,1.2)*exp(-9610/t);
}
double k50(double t, double sub_vec_nHx[]) {
	return 2.61e-10*exp(-8156/t);
}
double k51(double t, double sub_vec_nHx[]) {
	return 3.16e-10*exp(-21890/t);
}
double k52(double t, double sub_vec_nHx[]) {
	if (t <= 295) return 4.7e-11*pow(t/300,-0.34);
	else return 2.48e-12*pow(t/300,1.54)*exp(613/t);
}
double k53(double t, double sub_vec_nHx[]) {
	return 1.1e-10*sqrt(t/300)*exp(-77700/t);
}
double k54(double t, double sub_vec_nHx[]) {
	return 2.24e-9*pow(t/300,0.042)*exp(-t/466000);
}
double k55(double t, double sub_vec_nHx[]) {
	return 7.7e-9*exp(-17560/t);
}
double k56(double t, double sub_vec_nHx[]) {
	return 2.4e-9;
}
double k57(double t, double sub_vec_nHx[]) {
	return 2.0e-9;
}
double k58(double t, double sub_vec_nHx[]) {
	return 1e-10*exp(-4640/t);
}
double k59(double t, double sub_vec_nHx[]) {
	return 7.5e-10;
}
double k60(double t, double sub_vec_nHx[]) {
	return 1.2e-9;
}
double k61(double t, double sub_vec_nHx[]) {
	return 3.5e-10;
}
double k62(double t, double sub_vec_nHx[]) {
	return 1.4e-9;
}
double k63(double t, double sub_vec_nHx[]) {
	return 1e-9*exp(-7080/t);
}
double k64(double t, double sub_vec_nHx[]) {
	return 1.6e-9;
}
double k65(double t, double sub_vec_nHx[]) {
	return 7.5e-10;
}
double k66(double t, double sub_vec_nHx[]) {
	return 7.0e-10*exp(-10560/t);
}
double k67(double t, double sub_vec_nHx[]) {
	return 4.0e-10;
}
double k68(double t, double sub_vec_nHx[]) {
	return 4.8e-10;
}
double k69(double t, double sub_vec_nHx[]) {
	return 1.7e-9;
}
double k70(double t, double sub_vec_nHx[]) {
	return 1.5e-9;
}
double k71(double t, double sub_vec_nHx[]) {
	return 8.4e-10;
}
double k72(double t, double sub_vec_nHx[]) {
	return 1.3e-9;
}
double k73(double t, double sub_vec_nHx[]) {
	return 7.7e-10;
}
double k74(double t, double sub_vec_nHx[]) {
	return 1.01e-9;
}
double k75(double t, double sub_vec_nHx[]) {
	return 6.4e-10;
}
double k76(double t, double sub_vec_nHx[]) {
	return 5.9e-9;
}
double k77(double t, double sub_vec_nHx[]) {
	return 9.0e-10;
}
double k78(double t, double sub_vec_nHx[]) {
	return 1.8e-9;
}
double k79(double t, double sub_vec_nHx[]) {
	return 1e-11;
}
double k80(double t, double sub_vec_nHx[]) {
	return 3.8e-10;
}
double k81(double t, double sub_vec_nHx[]) {
	return 6.2e-10;
}
double k82(double t, double sub_vec_nHx[]) {
	return 9.1e-10;
}
double k83(double t, double sub_vec_nHx[]) {
	return 5.2e-11;
}
double k84(double t, double sub_vec_nHx[]) {
	return 2.7e-11;
}
double k85(double t, double sub_vec_nHx[]) {
	return 1.7e-9;
}
double k86(double t, double sub_vec_nHx[]) {
	return 1.1e-9;
}
double k87(double t, double sub_vec_nHx[]) {
	return 2.5e-9;
}
double k88(double t, double sub_vec_nHx[]) {
	return 7.2e-15;
}
double k89(double t, double sub_vec_nHx[]) {
	return 3.7e-14*exp(-35/t);
}
double k90(double t, double sub_vec_nHx[]) {
	return 1.9e-9;
}
double k91(double t, double sub_vec_nHx[]) {
	return 1.4e-9;
}
double k92(double t, double sub_vec_nHx[]) {
	return 7.5e-10;
}
double k93(double t, double sub_vec_nHx[]) {
	return 1.6e-9;
}
double k94(double t, double sub_vec_nHx[]) {
	return 2.1e-9;
}
double k95(double t, double sub_vec_nHx[]) {
	return 1.1e-9;
}
double k96(double t, double sub_vec_nHx[]) {
	return 6.9e-9;
}
double k97(double t, double sub_vec_nHx[]) {
	return 2.04e-10;
}
double k98(double t, double sub_vec_nHx[]) {
	return 2.86e-10;
}
double k99(double t, double sub_vec_nHx[]) {
	return 6.05e-11;
}
double k100(double t, double sub_vec_nHx[]) {
	return 2e-9;
}
double k101(double t, double sub_vec_nHx[]) {
	return 3.3e-11;
}
double k102(double t, double sub_vec_nHx[]) {
	return 1.1e-9;
}
double k103(double t, double sub_vec_nHx[]) {
	return 5.2e-11;
}
double k104(double t, double sub_vec_nHx[]) {
	return 1.4e-9*pow(t/300,-0.5);
}
double k105(double t, double sub_vec_nHx[]) {
	return 1.4e-16*pow(t/300,-0.5);
}
double k106(double t, double sub_vec_nHx[]) {
	return 7.5e-10;
}
double k107(double t, double sub_vec_nHx[]) {
	return 2.3e-7*pow(t/300,-0.5);
}
double k108(double t, double sub_vec_nHx[]) {
	return 2.3e-7*pow(t/300,-0.5);
}
double k109(double t, double sub_vec_nHx[]) {
	return 2.32e-7*pow(t/300,-0.52)*exp(t/22400);
}
double k110(double t, double sub_vec_nHx[]) {
	return 2.34e-8*pow(t/300,-0.52);
}
double k111(double t, double sub_vec_nHx[]) {
	return 4.36e-8*pow(t/300,-0.52);
}
double k112(double t, double sub_vec_nHx[]) {
	return 7e-8*pow(t/300,-0.5);
}
double k113(double t, double sub_vec_nHx[]) {
	return 1.6e-7*pow(t/300,-0.6);
}
double k114(double t, double sub_vec_nHx[]) {
	return 4.03e-7*pow(t/300,-0.6);
}
double k115(double t, double sub_vec_nHx[]) {
	return 7.68e-8*pow(t/300,-0.6);
}
double k116(double t, double sub_vec_nHx[]) {
	return 7.75e-8*pow(t/300,-0.5);
}
double k117(double t, double sub_vec_nHx[]) {
	return 1.95e-7*pow(t/300,-0.5);
}
double k118(double t, double sub_vec_nHx[]) {
	return 2.0e-7*pow(t/300,-0.4);
}
double k119(double t, double sub_vec_nHx[]) {
	return 6.3e-9*pow(t/300,-0.48);
}
double k120(double t, double sub_vec_nHx[]) {
	return 3.05e-7*pow(t/300,-0.5);
}
double k121(double t, double sub_vec_nHx[]) {
	return 3.9e-8*pow(t/300,-0.5);
}
double k122(double t, double sub_vec_nHx[]) {
	return 8.6e-8*pow(t/300,-0.5);
}
double k123(double t, double sub_vec_nHx[]) {
	return 1.08e-7*pow(t/300,-0.5);
}
double k124(double t, double sub_vec_nHx[]) {
	return 6.02e-8*pow(t/300,-0.5);
}
double k125(double t, double sub_vec_nHx[]) {
	return 2.58e-7*pow(t/300,-0.5);
}
double k126(double t, double sub_vec_nHx[]) {
	return 5.6e-9*pow(t/300,-0.5);
}
double k127(double t, double sub_vec_nHx[]) {
	return 1.95e-7*pow(t/300,-0.7);
}
double k128(double t, double sub_vec_nHx[]) {
	return 2.75e-7*pow(t/300,-0.55);
}
double k129(double t, double sub_vec_nHx[]) {
	return 2.76e-7*pow(t/300,-0.64);
}
double k130(double t, double sub_vec_nHx[]) {
	return 2.4e-8*pow(t/300,-0.64);
}
double k131(double t, double sub_vec_nHx[]) {
	return 1.1e-7*pow(t/300,-1.0);
}
double k132(double t, double sub_vec_nHx[]) {
	return 1e-9;
}
double k133(double t, double sub_vec_nHx[]) {
	return 1e-9;
}
double k134(double t, double sub_vec_nHx[]) {
	return 1e-10;
}
double k135(double t, double sub_vec_nHx[]) {
	return 5e-10;
}
double k136(double t, double sub_vec_nHx[]) {
	return 1e-13;
}
double k137(double t, double sub_vec_nHx[]) {
	return 5e-10;
}
double k138(double t, double sub_vec_nHx[]) {
	return 5e-10;
}
double k139(double t, double sub_vec_nHx[]) {
	return 7e-10;
}
double k140(double t, double sub_vec_nHx[]) {
	return 5e-10;
}
double k141(double t, double sub_vec_nHx[]) {
	return 1e-16;
}
double k142(double t, double sub_vec_nHx[]) {
	return 2.25e-15;
}
double k143(double t, double sub_vec_nHx[]) {
	return 1e-17;
}
double k144(double t, double sub_vec_nHx[]) {
	return 1e-17;
}
double k145(double t, double sub_vec_nHx[]) {
	return 4.36e-18*pow(t/300,0.35)*exp(-161.3/t);
}
double k146(double t, double sub_vec_nHx[]) {
	if (t <= 300) return 2.1e-19;
	else return 3.09e-17*pow(t/300,0.33)*exp(-1629/t);
}
double k147(double t, double sub_vec_nHx[]) {
	return 4.46e-16*pow(t,-0.5)*exp(-4.93/pow(t,0.6667));
}
double k148(double t, double sub_vec_nHx[]) {
	return 4e-16*pow(t/300,-0.2);
}
double k149(double t, double sub_vec_nHx[]) {
	if (t <= 300) return  2.5e-18;
	else return 3.14e-18*pow(t/300,-0.15)*exp(68/t);
}
double k150(double t, double sub_vec_nHx[]) {
	return 1.5e-15;
}
double k151(double t, double sub_vec_nHx[]) {
	return 9.9e-19*pow(t/300,-0.38);
}
double k152(double t, double sub_vec_nHx[]) {
	return 4.9e-20*pow(t/300,1.58);
}
double k153(double t, double sub_vec_nHx[]) {
	return 5.26e-18*pow(t/300,-5.22)*exp(-90/t);
}
double k154(double t, double sub_vec_nHx[]) {
	if (t <= 300) return 1.32e-32*pow(t/300,-0.38);
	else return 1.32e-32*pow(t/300,-1.0);
}
double k155(double t, double sub_vec_nHx[]) {
	return 2.8e-31*pow(t,-0.6);
}
double k156(double t, double sub_vec_nHx[]) {
	return 6.9e-32*pow(t,-0.4);
}
double k157(double t, double sub_vec_nHx[]) {
	if (t <= 5000) return 5.99e-33*pow(t/5000,-1.6);
	else return 5.99e-33*pow(t/5000,-0.64)*exp(5255/t);
}
double k158(double t, double sub_vec_nHx[]) {
	if (t < 2000) return 6.16e-29*pow(t/300,-3.08);
	else return 2.14e-29*pow(t/300,-3.08)*exp(2114/t);
}
double k159(double t, double sub_vec_nHx[]) {
	// return 100*k158(t);
	return 1e-19*pow(t,-3.08)*exp(2114/t);
}
double k160(double t, double sub_vec_nHx[]) {
	// return 100*k158(t);
	return 1e-19*pow(t,-3.08)*exp(2114/t);
}
double k161(double t, double sub_vec_nHx[]) {
	return 4.33e-32*pow(t/300,-1.0);
}
double k162(double t, double sub_vec_nHx[]) {
	return 2.56e-31*pow(t/300,-2.0);
}
double k163(double t, double sub_vec_nHx[]) {
	return 9.2e-34*pow(t/300,-1.0);
}
double k164(double t, double sub_vec_nHx[]) {
	return 2.0e-11*pow(t/300,0.44);
}
double k165(double t, double temp_grain) {
	double fA = 1/(1.0 + 1e4*exp(-600/temp_grain));
	return (3.0e-18*sqrt(t))/(fA*(1.0 + 0.04*sqrt(t + temp_grain) + 0.002*t + 8e-6*t*t));
}
