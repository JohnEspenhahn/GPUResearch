package com.hahn;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;

public class Main {
	private static List<Equation> equations;
	private static Map<String, ODEEntry> ode_eqs;

	public static void main(String[] args) throws IOException {
		equations = new ArrayList<Equation>();
		ode_eqs = new HashMap<String, ODEEntry>();
		ode_eqs.put("H_p", new ODEEntry("H_p"));
		ode_eqs.put("H2", new ODEEntry("H2"));
		ode_eqs.put("He_p", new ODEEntry("He_p"));
		ode_eqs.put("C_p", new ODEEntry("C_p"));
		ode_eqs.put("O_p", new ODEEntry("O_p"));
		ode_eqs.put("OH", new ODEEntry("OH"));
		ode_eqs.put("H2O", new ODEEntry("H2O"));
		ode_eqs.put("CO", new ODEEntry("CO"));
		ode_eqs.put("C2", new ODEEntry("C2"));
		ode_eqs.put("O2", new ODEEntry("O2"));
		ode_eqs.put("HCO_p", new ODEEntry("HCO_p"));
		ode_eqs.put("CH", new ODEEntry("CH"));
		ode_eqs.put("CH2", new ODEEntry("CH2"));
		ode_eqs.put("CH3_p", new ODEEntry("CH3_p"));
		
		// Equilibrium
		ode_eqs.put("H_m", new ODEEntry("H_m"));
		ode_eqs.put("H2_p", new ODEEntry("H2_p"));
		ode_eqs.put("H3_p", new ODEEntry("H3_p"));
		ode_eqs.put("CH_p", new ODEEntry("CH_p"));
		ode_eqs.put("CH2_p", new ODEEntry("CH2_p"));
		ode_eqs.put("OH_p", new ODEEntry("OH_p"));
		ode_eqs.put("H2O_p", new ODEEntry("H2O_p"));
		ode_eqs.put("H3O_p", new ODEEntry("H3O_p"));
		ode_eqs.put("CO_p", new ODEEntry("CO_p"));
		ode_eqs.put("HOC_p", new ODEEntry("HOC_p"));
		ode_eqs.put("O_m", new ODEEntry("O_m"));
		ode_eqs.put("C_m", new ODEEntry("C_m"));
		ode_eqs.put("O2_p", new ODEEntry("O2_p"));
		
		
		Scanner s = new Scanner(new File("equations.csv"));
		
		int i = 1;
		while (s.hasNextLine()) {
			String streq = s.nextLine();
			try {
				Equation eq = new Equation(i, streq);
				equations.add(eq);
				
				for (Entry<String, ODEEntry> ode: ode_eqs.entrySet()) {
					int contains = eq.contains(ode.getKey());
					if (contains > 0) {
						ode.getValue().addCreationEq(eq);
					} else if (contains < 0) {
						ode.getValue().addDistructionEq(eq);
					}
				}
				
				// System.out.println(eq);
			} catch (Exception e) {
				System.err.println(i + " " + streq);
			}
			
			i += 1;
		}
		
		s.close();
		
		ode_out();
		k_out();
	}
	
	private static void ode_out() throws IOException {
		File f = new File("ode_out.txt");
		if (f.exists()) f.delete();
		
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f.getCanonicalPath(), true)));
		for (Entry<String, ODEEntry> ode: ode_eqs.entrySet()) {
			out.println(ODEEntryFormat.formatShort(ode.getValue()));	    
		}
		out.close();
	}
	
	private static void k_out() throws IOException {
		File f = new File("k_out.txt");
		if (f.exists()) f.delete();
		
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f.getCanonicalPath(), true)));
		for (Equation e: equations) {
			out.println(EquationFormat.format(e));
		}
		out.close();
	}
	
}
