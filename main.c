#include "simulation.h"

char *p, s[100];

double getd(double, bool);

int main() {	
	do {
		printf("Number of seconds to run (0 to quit): ");
		double x2 = getd(0, false);
		
		if (x2 > 0) {
			printf("Accuracy eps (default 1): ");
			int eps = getd(1, true);
			
			printf("Initial setp (default 100): ");
			int h1 = getd(100, true);
			
			FILE *fp = fopen("out.csv", "w");
			
			run(fp, x2, eps, h1, true);
			
			fclose(fp);
		} else break;
	} while (1);
	
	return 0;
}

double getd(double def, bool has_def) {
	double i = 0;
	p = 0;
	while (fgets(s, sizeof(s), stdin)) {
		i = strtod(s, &p);
		if (p == s || *p != '\n') {
			if (has_def) return def;
			printf("Please enter an integer: ");
		} else break;
	}
	return i;
}