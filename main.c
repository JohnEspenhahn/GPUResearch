#include "simulation.h"

int main() {
	FILE *fp = fopen("out.csv", "w");
	
	char *p, s[100];
	double x2;
	
	do {
		printf("Number of seconds to run (0 to quit): ");
		p = 0;

		while (fgets(s, sizeof(s), stdin)) {
			x2 = strtod(s, &p);
			if (p == s || *p != '\n') {
				printf("Please enter an integer: ");
			} else break;
		}
		
		if (x2 > 0) run(fp, x2, true);
		else break;
	} while (1);
	
	fclose(fp);
	
	return 0;
}