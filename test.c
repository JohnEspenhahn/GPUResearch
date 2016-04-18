#include "test.h"
#include "g10/derivs.h"
#include "odeint.h"
#include "stiff.h"
#include <stdio.h>
#include "simulation.h"
#include "jacobian.h"

double end_pH2 = 1e-22;
double end_pH_p = 1e-25;
double end_temp = 199;

double (*karr[])(double) = {
k1,
k2,
k3,
k4,
k5,
k6,
k7,
k8,
// k9,
// k10,
k11,
k12,
k13,
k14,
k15,
k16,
k17,
k18,
k19,
k20,
k21,
k22,
k23,
k24,
k25,
k26,
k27,
k28,
k29,
// k30,
k31,
k32,
k33,
k34,
k35,
k36,
k37,
k38,
k39,
k40,
k41,
k42,
k43,
k44,
k45,
k46,
k47,
k48,
k49,
k50,
k51,
k52,
k53,
k54,
k55,
k56,
k57,
k58,
k59,
k60,
k61,
k62,
k63,
k64,
k65,
k66,
k67,
k68,
k69,
k70,
k71,
k72,
k73,
k74,
k75,
k76,
k77,
k78,
k79,
k80,
k81,
k82,
k83,
k84,
k85,
k86,
k87,
k88,
k89,
k90,
k91,
k92,
k93,
k94,
k95,
k96,
k97,
k98,
k99,
k100,
k101,
k102,
k103,
k104,
k105,
k106,
k107,
k108,
k109,
k110,
k111,
k112,
k113,
k114,
k115,
k116,
k117,
k118,
k119,
k120,
k121,
k122,
k123,
k124,
k125,
k126,
k127,
k128,
k129,
k130,
k131,
k132,
k133,
k134,
k135,
k136,
k137,
k138,
k139,
k140,
k141,
k142,
k143,
k144,
k145,
k146,
k147,
k148,
k149,
k150,
k151,
k152,
k153,
k154,
k155,
k156,
k157,
k158,
k159,
k160,
k161,
k162,
k163,
k164,
// k165
};

void testks(FILE *fp, double t) {
	for (int i = 0; i < sizeof(karr)/sizeof(double*); i++) {
		fprintf(fp, "%d,%G\n", i, karr[i](t));
	}
}

void outputRates() {
	FILE *fp = fopen("out.csv", "w");
	
	int steps = 5;	
	
	testks(fp, 50.0);
	
	// compareJacobn(fp, steps);
	// testCopy();
	
	fclose(fp);
}

int main() {
	outputRates();
	
	return 0;
}