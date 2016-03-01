package com.hahn;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Equation {
	private static final Pattern elements = Pattern.compile("([A-Za-z][a-z]*)([0-9]?[\\-\\+]?)");
	
	int id;
	String streq;
	Map<String, Byte> reactants;
	Map<String, Byte> products;
	Map<String, Byte> net_products;
	
	public Equation(int id, String streq) {
		this.id = id;
		this.streq = streq;
		this.reactants = new HashMap<String, Byte>();
		this.products = new HashMap<String, Byte>();
		this.net_products = new HashMap<String, Byte>();
		
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
	
	private void parseElements(String eq, Map<String, Byte> map, int species_add, Map<String, Byte> species_map) {
		for (String entry: eq.split("\\s+")) {
			Matcher matcher = elements.matcher(entry);
			while (matcher.find()) {
				String species = matcher.group();
				String element = matcher.group(1);
				String count   = matcher.group(2);
				
				// Track species
				Byte oldSpeciesCount = species_map.get(species);
				if (oldSpeciesCount == null) oldSpeciesCount = 0;
				species_map.put(species, (byte) (oldSpeciesCount + species_add));
				
				// Check charge
				if (count != null && (count.endsWith("-") || count.endsWith("+"))) {
					Byte oldCharge = map.get("e");
					if (oldCharge == null) oldCharge = 0;
					
					if (count.endsWith("-")) map.put("e", (byte) (oldCharge + 1));
					else map.put("e", (byte) (oldCharge - 1));
					
					count = count.substring(0, count.length()-1);
				}
				
				// Check count
				int addCount = 1;
				if (count != null && count.length() > 0) {
					addCount = Integer.parseInt(count);
				}
				
				// Track element
				Byte oldCount = map.get(element);
				if (oldCount == null) oldCount = 0;
				map.put(element, (byte) (oldCount + addCount));
			}
		}
	}
	
	public String toString() {
		String str = "";
		for (Entry<String, Byte> element: net_products.entrySet()) {
			str += element.getKey() + ":" + element.getValue() + ",";
		}
		
		return id + " " + streq + " ; " + str;
	}
}
