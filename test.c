#include "main.h"
#include <stdio.h>

int main() {
	printf("k2(600, 91, 0): %G\n", k2(600, 91, 0));
	printf("k2_old(600, 91, 0): %G\n", k2_old(500, 91, 0));
	printf("k3(600, 91, 0): %G\n", k3(500, 91, 0));
	
	return 0;
}