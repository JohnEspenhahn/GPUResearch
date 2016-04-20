package com.hahn;

public class ODEEntryFormat {

	public static String format(ODEEntry e) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("double dn%s(double sub_vec_nHx[]) {\n", e.getDeriv()));
		sb.append("	double nH_p = sub_vec_nHx[0], nH2  = sub_vec_nHx[1], nHe_p = sub_vec_nHx[2]\n");
		sb.append("			 , nC_p = sub_vec_nHx[3], nO_p = sub_vec_nHx[4], nOH   = sub_vec_nHx[5]\n");
		sb.append("			 , nH2O = sub_vec_nHx[6], nCO  = sub_vec_nHx[7], nC2   = sub_vec_nHx[8]\n");
		sb.append("			 , nO2  = sub_vec_nHx[9], nHCO_p = sub_vec_nHx[10], nCH = sub_vec_nHx[11]\n");
		sb.append("			 , nCH2 = sub_vec_nHx[12], nCH3_p = sub_vec_nHx[13], nH_m = sub_vec_nHx[14]\n");
		sb.append("			 , nH2_p = sub_vec_nHx[15], nH3_p = sub_vec_nHx[16], nCH_p = sub_vec_nHx[17]\n");
		sb.append("			 , nCH2_p = sub_vec_nHx[18], nOH_p = sub_vec_nHx[19], nH2O_p = sub_vec_nHx[20]\n");
		sb.append("			 , nH3O_p = sub_vec_nHx[21], nCO_p = sub_vec_nHx[22], nHOC_p = sub_vec_nHx[23]\n");
		sb.append("			 , nO_m = sub_vec_nHx[24], nC_m = sub_vec_nHx[25], nO2_p = sub_vec_nHx[26]\n");
		sb.append("			 , nH = getnH(nH_p, nH2, nOH, nH2O, nHCO_p, nCH, nCH2, nCH3_p)\n");
		sb.append("			 , nC = getnC(nC_p, nCO, nC2, nHCO_p, nCH, nCH2, nCH3_p)\n");
		sb.append("			 , ne = getne(nH_p, nHe_p, nC_p, nO_p, nHCO_p, nCH3_p)\n");
		sb.append("			 , nO = getnO(nO_p, nOH, nH2O, nCO, nO2, nHCO_p)\n");
		sb.append("			 , nHe = getnHe(nHe_p)\n");
		sb.append("			 , nM = nC + nO + nSi\n");
		sb.append("			 , t = T;\n\n");
		sb.append(String.format("\treturn %s;\n", e.toString()));
		sb.append("}\n");
		
		return sb.toString();
	}
	
}
