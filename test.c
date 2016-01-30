#include "main.h"
#include <stdio.h>

#define TEMP 1000 // K
#define GRAIN_TEMP 75 // K

int main() {
	printf("k1(%d, %d): %G\n", TEMP, GRAIN_TEMP, k1(TEMP, GRAIN_TEMP));
	printf("k2(%d, %d): %G\n", TEMP, GRAIN_TEMP, -k2(TEMP, 91, 0));
	printf("k3(%d, %d): %G\n", TEMP, GRAIN_TEMP, -k3(TEMP, 91, 0));
	
	printf("k6(%d, %d): %G\n", TEMP, GRAIN_TEMP, k6(TEMP));
	printf("k7(%d, %d): %G\n", TEMP, GRAIN_TEMP, -k7(TEMP));
	printf("k8(%d, %d): %G\n", TEMP, GRAIN_TEMP, -k8(TEMP));
	
	return 0;
}