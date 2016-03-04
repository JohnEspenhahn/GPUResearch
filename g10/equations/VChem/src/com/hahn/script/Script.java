package com.hahn.script;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Script {

	public static void main(String[] args) throws IOException {
		File f = new File("k.txt");
		if (f.exists()) f.delete();
		
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f.getCanonicalPath(), true)));
		for (int i = 1; i <= 165; i++) {
			out.println("double k" + i + "(double t) {");
			out.println("\treturn ");
			out.println("}");
		}
		out.close();
	}

}
