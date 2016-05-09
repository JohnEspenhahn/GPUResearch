package com.hahn;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Equation {
	private static final Pattern elements = Pattern.compile("([A-Za-z][a-z]*)([0-9]?_?[mp]?)");
	
	private int id;
	private String streq;
	private Map<String, Byte> reactants;
	private Map<String, Byte> products;
	private Map<String, EquationProduct> net_products;
	
	public Equation(int id, String streq) {
		this.id = id;
		this.streq = streq;
		this.reactants = new HashMap<String, Byte>();
		this.products = new HashMap<String, Byte>();
		this.net_products = new HashMap<String, EquationProduct>();
		
		String[] rp = streq.trim().split("=");
		parseElements(rp[0], reactants, -1, net_products);
		parseElements(rp[1], products ,  1, net_products);
		validateElements(reactants, products);
		validateElements(products, reactants);
	}
	
	private void validateElements(Map<String, Byte> a, Map<String, Byte> b) { 
		for (Entry<String, Byte> element: a.entrySet()) {
			if (element.getKey().equals("g")) continue;
			
			Byte cnt = b.get(element.getKey());
			if (cnt == null) cnt = 0;
			
			if (cnt != element.getValue()) {
				throw new RuntimeException(String.format("Unballanced equation '%s' for element '%s'", streq, element.getKey()));
			}
		}
	}
	
	private void parseElements(String eq, Map<String, Byte> map, int species_add, Map<String, EquationProduct> species_map) {
		for (String species: eq.split("\\s+")) {
			if (species.length() == 0) continue;
			
			// Track species
			EquationProduct speciesCount = species_map.getOrDefault(species, new EquationProduct(species));
			speciesCount.add(species_add);
			species_map.put(species, speciesCount);
			
			// Track elements
			Matcher matcher = elements.matcher(species);
			while (matcher.find()) {
				String element = matcher.group(1);
				String count   = matcher.group(2);
				
				// Check charge
				if (count != null && (count.endsWith("m") || count.endsWith("p"))) {
					Byte oldCharge = map.getOrDefault("e", (byte) 0);
					
					if (count.endsWith("m")) map.put("e", (byte) (oldCharge + 1));
					else map.put("e", (byte) (oldCharge - 1));
					
					count = count.substring(0, count.length()-2);
				}
				
				// Check count
				int addCount = 1;
				if (count != null && count.length() > 0) {
					addCount = Integer.parseInt(count);
				}
				
				// Track element
				Byte oldCount = map.getOrDefault(element, (byte) 0);
				map.put(element, (byte) (oldCount + addCount));
			}
		}
	}
	
	public int getId() {
		return id;
	}
	
	public int contains(String species) {
		EquationProduct p = net_products.get(species);
		return (p == null ? 0 : p.getNetProduction());
	}
	
	public String asODEEntry() {
		String str = "k" + id + "(t)";
		for (Entry<String, EquationProduct> e: net_products.entrySet()) {
			EquationProduct p = e.getValue();
			if (p.isReactant()) {
				for (int j = 0; j < p.getDistruction(); j++) {
					str += "*n" + e.getKey();
				}
			}
		}
		
		return str;
	}
	
	/**
	 * Get the base-1 indexed reactant
	 * @param idx The base-1 index
	 * @return The reaction
	 * @throws IndexOutOfBoundsException If no idx'th reaction
	 */
	public EquationProduct getNReactant(int idx) {
		int at_idx = 1;
		for (Entry<String, EquationProduct> e: net_products.entrySet()) {
			EquationProduct p = e.getValue();
			if (p.isReactant()) {
				for (int j = 0; j < p.getDistruction(); j++) {
					at_idx += 1;
					if (at_idx > idx) return p;
				}
			}
		}
		
		throw new IndexOutOfBoundsException("No " + at_idx + "th reactant in reaction " + getId());
	}
	
	public String toString() {
		String str = "";
		for (Entry<String, EquationProduct> element: net_products.entrySet()) {
			str += element.getKey() + ":" + element.getValue().getNetProduction() + ",";
		}
		
		return "" + id + " " + streq + " ; " + str;
	}
}
