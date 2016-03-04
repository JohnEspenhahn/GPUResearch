package com.hahn;

import java.util.ArrayList;
import java.util.List;

public class ODEEntry {
	private List<Equation> c, d;
	
	public ODEEntry() { 
		c = new ArrayList<Equation>();
		d = new ArrayList<Equation>();
	}
	
	public void addCreationEq(Equation c) {
		this.c.add(c);
	}
	
	public void addDistructionEq(Equation d) {
		this.d.add(d);
	}
	
	@Override
	public String toString() {
		String str = "";
		boolean first = true;
		for (Equation e: c) {
			if (first) first = false;
			else str += " + ";
			str += e.asODEEntry();
		}
		
		for (Equation e: d) {
			str += " - " + e.asODEEntry();
		}
		
		return str;
	}
}
