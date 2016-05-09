package com.hahn;

public class EquationProduct {
	private int creation, distruction;
	private String species;
	
	public EquationProduct(String species) {
		this.species = species;
	}
	
	public String getSpecies() {
		return this.species;
	}
	
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
	
	@Override
	public String toString() {
		 return getSpecies(); // String.format("C:%d,D:%d", getCreation(), getDistruction());
	}
}
