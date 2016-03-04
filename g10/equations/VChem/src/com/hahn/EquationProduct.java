package com.hahn;

public class EquationProduct {
	private int creation, distruction;
	
	public void add(int amnt) {
		if (amnt > 0) {
			creation += amnt;
		} else if (amnt < 0) {
			distruction += -amnt;
		}
	}
	
	public int getCreation() {
		return creation;
	}
	
	public int getDistruction() {
		return distruction;
	}
	
	public int getNetProduction() {
		return creation - distruction;
	}
	
	public boolean isReactant() {
		return getDistruction() > 0;
	}
}
