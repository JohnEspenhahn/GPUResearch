package com.hahn;

public class EquationFormat {

	public static String format(Equation e) {
		String str = String.format("{ .k = &k%d, .n1 = %s, .n2 = %s },", e.getId(), encode(e.getNReactant(1)), encode(e.getNReactant(2)));
		
		try {
			e.getNReactant(3);
			System.err.println("WARNING: equation " + e.getId() + " has 3 or more nreactants!");
		} catch (IndexOutOfBoundsException exc) {}
		
		return str;
	}
	
	private static String encode(EquationProduct p) {
		switch (p.getSpecies()) {
		case "H":  return "100";
		case "He": return "101";
		case "C":  return "102";
		case "O":  return "103";
		case "e":  return "104";
		case "M":  return "105";
		default: return "n" + p.getSpecies();
		}
	}
}
