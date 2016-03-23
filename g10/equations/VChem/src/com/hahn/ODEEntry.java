package com.hahn;

import java.util.ArrayList;
import java.util.List;

public class ODEEntry {
	private List<Equation> creation, distruction;
	private String derv;
	
	public ODEEntry(String derv) {
		this.derv = derv;
		this.creation = new ArrayList<Equation>();
		this.distruction = new ArrayList<Equation>();
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
			str += e.asODEEntry(derv.replace("_p", ""));
		}
		
		for (Equation e: distruction) {
			str += " - " + e.asODEEntry(derv);
		}
		
		return str;
	}
}
