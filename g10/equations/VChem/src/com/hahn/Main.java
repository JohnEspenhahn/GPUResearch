package com.hahn;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

public class Main {
	static List<Equation> equations;
	static Map<String, List<Equation>> ode_eqs;

	public static void main(String[] args) throws FileNotFoundException {
		equations = new ArrayList<Equation>();
		ode_eqs = new HashMap<String, List<Equation>>();
		
		File f = new File("equations.csv");
		Scanner s = new Scanner(f);
		
		int i = 1;
		while (s.hasNextLine()) {
			String streq = s.nextLine();
			try {
				Equation eq = new Equation(i, streq);
				equations.add(eq);
				
				System.out.println(eq);
			} catch (Exception e) {
				System.err.println(i + " " + streq);
			}
			
			i += 1;
		}
		
		s.close();
	}
	
}
