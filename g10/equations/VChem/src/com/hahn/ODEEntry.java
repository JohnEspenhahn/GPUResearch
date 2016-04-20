package com.hahn;

import java.util.ArrayList;
import java.util.List;

public class ODEEntry {
	private List<Equation> creation, distruction;
	
	/** The species for which this is determining the derivative */
	private String derv;
	
	public ODEEntry(String derv) {
		this.derv = derv;
		this.creation = new ArrayList<Equation>();
		this.distruction = new ArrayList<Equation>();
	}
	
	public String getDeriv() {
		return derv;
	}
	
	public void addCreationEq(Equation c) {
		this.creation.add(c);
	}
	
	public void addDistructionEq(Equation d) {
		this.distruction.add(d);
	}
	
	@Override
	public String toString() {
		String str = "";
		boolean first = true;
		for (Equation e: creation) {
			if (first) first = false;
			else str += " + ";
			str += e.asODEEntry();
		}
		
		for (Equation e: distruction) {
			str += " - " + e.asODEEntry();
		}
		
		return str;
	}
}
