#include "tscalecsv.h"

struct reaction reactions[] = {
{ .k = &k1, .n1 = 104, .n2 = 100 },
{ .k = &k2, .n1 = 15, .n2 = 100 },
{ .k = &k3, .n1 = 100, .n2 = 1 },
{ .k = &k4, .n1 = 100, .n2 = 16 },
{ .k = &k5, .n1 = 15, .n2 = 1 },
{ .k = &k6, .n1 = 104, .n2 = 16 },
{ .k = &k7, .n1 = 1, .n2 = 2 },
{ .k = &k8, .n1 = 104, .n2 = 2 },
{ .k = &k9, .n1 = 100, .n2 = 2 },
{ .k = &k10, .n1 = 2, .n2 = 2 },
{ .k = &k11, .n1 = 104, .n2 = 100 },
{ .k = &k12, .n1 = 104, .n2 = 1 },
{ .k = &k13, .n1 = 104, .n2 = 15 },
{ .k = &k14, .n1 = 15, .n2 = 100 },
{ .k = &k15, .n1 = 15, .n2 = 1 },
{ .k = &k16, .n1 = 104, .n2 = 101 },
{ .k = &k17, .n1 = 104, .n2 = 3 },
{ .k = &k18, .n1 = 100, .n2 = 3 },
{ .k = &k19, .n1 = 1, .n2 = 101 },
{ .k = &k20, .n1 = 4, .n2 = 104 },
{ .k = &k21, .n1 = 5, .n2 = 104 },
{ .k = &k22, .n1 = 102, .n2 = 104 },
{ .k = &k23, .n1 = 104, .n2 = 103 },
{ .k = &k24, .n1 = 5, .n2 = 100 },
{ .k = &k25, .n1 = 1, .n2 = 103 },
{ .k = &k26, .n1 = 3, .n2 = 103 },
{ .k = &k27, .n1 = 102, .n2 = 1 },
{ .k = &k28, .n1 = 4, .n2 = 100 },
{ .k = &k29, .n1 = 102, .n2 = 3 },
{ .k = &k30, .n1 = 2, .n2 = 101 },
{ .k = &k31, .n1 = 100, .n2 = 6 },
{ .k = &k32, .n1 = 24, .n2 = 2 },
{ .k = &k33, .n1 = 24, .n2 = 8 },
{ .k = &k34, .n1 = 102, .n2 = 2 },
{ .k = &k35, .n1 = 12, .n2 = 100 },
{ .k = &k36, .n1 = 12, .n2 = 2 },
{ .k = &k37, .n1 = 102, .n2 = 12 },
{ .k = &k38, .n1 = 12, .n2 = 103 },
{ .k = &k39, .n1 = 100, .n2 = 13 },
{ .k = &k40, .n1 = 13, .n2 = 103 },
{ .k = &k41, .n1 = 13, .n2 = 103 },
{ .k = &k42, .n1 = 9, .n2 = 103 },
{ .k = &k43, .n1 = 2, .n2 = 103 },
{ .k = &k44, .n1 = 100, .n2 = 6 },
{ .k = &k45, .n1 = 6, .n2 = 2 },
{ .k = &k46, .n1 = 102, .n2 = 6 },
{ .k = &k47, .n1 = 6, .n2 = 103 },
{ .k = &k48, .n1 = 6, .n2 = 6 },
{ .k = &k49, .n1 = 7, .n2 = 100 },
{ .k = &k50, .n1 = 10, .n2 = 100 },
{ .k = &k51, .n1 = 10, .n2 = 2 },
{ .k = &k52, .n1 = 10, .n2 = 102 },
{ .k = &k53, .n1 = 100, .n2 = 8 },
{ .k = &k54, .n1 = 16, .n2 = 2 },
{ .k = &k55, .n1 = 100, .n2 = 17 },
{ .k = &k56, .n1 = 102, .n2 = 16 },
{ .k = &k57, .n1 = 102, .n2 = 17 },
{ .k = &k58, .n1 = 4, .n2 = 2 },
{ .k = &k59, .n1 = 100, .n2 = 18 },
{ .k = &k60, .n1 = 18, .n2 = 2 },
{ .k = &k61, .n1 = 18, .n2 = 103 },
{ .k = &k62, .n1 = 1, .n2 = 13 },
{ .k = &k63, .n1 = 100, .n2 = 19 },
{ .k = &k64, .n1 = 2, .n2 = 19 },
{ .k = &k65, .n1 = 19, .n2 = 103 },
{ .k = &k66, .n1 = 100, .n2 = 14 },
{ .k = &k67, .n1 = 14, .n2 = 103 },
{ .k = &k68, .n1 = 5, .n2 = 9 },
{ .k = &k69, .n1 = 5, .n2 = 2 },
{ .k = &k70, .n1 = 16, .n2 = 103 },
{ .k = &k71, .n1 = 17, .n2 = 103 },
{ .k = &k72, .n1 = 6, .n2 = 17 },
{ .k = &k73, .n1 = 4, .n2 = 6 },
{ .k = &k74, .n1 = 2, .n2 = 20 },
{ .k = &k75, .n1 = 21, .n2 = 2 },
{ .k = &k76, .n1 = 7, .n2 = 17 },
{ .k = &k77, .n1 = 7, .n2 = 4 },
{ .k = &k78, .n1 = 7, .n2 = 4 },
{ .k = &k79, .n1 = 102, .n2 = 22 },
{ .k = &k80, .n1 = 10, .n2 = 4 },
{ .k = &k81, .n1 = 10, .n2 = 4 },
{ .k = &k82, .n1 = 10, .n2 = 19 },
{ .k = &k83, .n1 = 27, .n2 = 102 },
{ .k = &k84, .n1 = 8, .n2 = 17 },
{ .k = &k85, .n1 = 8, .n2 = 17 },
{ .k = &k86, .n1 = 11, .n2 = 102 },
{ .k = &k87, .n1 = 11, .n2 = 7 },
{ .k = &k88, .n1 = 2, .n2 = 3 },
{ .k = &k89, .n1 = 2, .n2 = 3 },
{ .k = &k90, .n1 = 12, .n2 = 1 },
{ .k = &k91, .n1 = 1, .n2 = 13 },
{ .k = &k92, .n1 = 13, .n2 = 3 },
{ .k = &k93, .n1 = 3, .n2 = 9 },
{ .k = &k94, .n1 = 1, .n2 = 6 },
{ .k = &k95, .n1 = 6, .n2 = 3 },
{ .k = &k96, .n1 = 7, .n2 = 1 },
{ .k = &k97, .n1 = 7, .n2 = 3 },
{ .k = &k98, .n1 = 7, .n2 = 3 },
{ .k = &k99, .n1 = 7, .n2 = 3 },
{ .k = &k100, .n1 = 10, .n2 = 1 },
{ .k = &k101, .n1 = 10, .n2 = 3 },
{ .k = &k102, .n1 = 10, .n2 = 3 },
{ .k = &k103, .n1 = 27, .n2 = 102 },
{ .k = &k104, .n1 = 8, .n2 = 3 },
{ .k = &k105, .n1 = 8, .n2 = 3 },
{ .k = &k106, .n1 = 23, .n2 = 100 },
{ .k = &k107, .n1 = 26, .n2 = 1 },
{ .k = &k108, .n1 = 1, .n2 = 25 },
{ .k = &k109, .n1 = 15, .n2 = 3 },
{ .k = &k110, .n1 = 104, .n2 = 17 },
{ .k = &k111, .n1 = 104, .n2 = 17 },
{ .k = &k112, .n1 = 104, .n2 = 18 },
{ .k = &k113, .n1 = 104, .n2 = 19 },
{ .k = &k114, .n1 = 104, .n2 = 19 },
{ .k = &k115, .n1 = 104, .n2 = 19 },
{ .k = &k116, .n1 = 104, .n2 = 14 },
{ .k = &k117, .n1 = 104, .n2 = 14 },
{ .k = &k118, .n1 = 104, .n2 = 14 },
{ .k = &k119, .n1 = 104, .n2 = 20 },
{ .k = &k120, .n1 = 104, .n2 = 21 },
{ .k = &k121, .n1 = 104, .n2 = 21 },
{ .k = &k122, .n1 = 104, .n2 = 21 },
{ .k = &k123, .n1 = 104, .n2 = 22 },
{ .k = &k124, .n1 = 104, .n2 = 22 },
{ .k = &k125, .n1 = 104, .n2 = 22 },
{ .k = &k126, .n1 = 104, .n2 = 22 },
{ .k = &k127, .n1 = 27, .n2 = 104 },
{ .k = &k128, .n1 = 23, .n2 = 104 },
{ .k = &k129, .n1 = 11, .n2 = 104 },
{ .k = &k130, .n1 = 11, .n2 = 104 },
{ .k = &k131, .n1 = 104, .n2 = 24 },
{ .k = &k132, .n1 = 102, .n2 = 15 },
{ .k = &k133, .n1 = 15, .n2 = 103 },
{ .k = &k134, .n1 = 15, .n2 = 6 },
{ .k = &k135, .n1 = 26, .n2 = 100 },
{ .k = &k136, .n1 = 26, .n2 = 2 },
{ .k = &k137, .n1 = 26, .n2 = 103 },
{ .k = &k138, .n1 = 100, .n2 = 25 },
{ .k = &k139, .n1 = 2, .n2 = 25 },
{ .k = &k140, .n1 = 102, .n2 = 25 },
{ .k = &k141, .n1 = 1, .n2 = 2 },
{ .k = &k142, .n1 = 102, .n2 = 104 },
{ .k = &k143, .n1 = 102, .n2 = 100 },
{ .k = &k144, .n1 = 102, .n2 = 2 },
{ .k = &k145, .n1 = 102, .n2 = 102 },
{ .k = &k146, .n1 = 102, .n2 = 103 },
{ .k = &k147, .n1 = 4, .n2 = 100 },
{ .k = &k148, .n1 = 4, .n2 = 2 },
{ .k = &k149, .n1 = 4, .n2 = 103 },
{ .k = &k150, .n1 = 104, .n2 = 103 },
{ .k = &k151, .n1 = 100, .n2 = 103 },
{ .k = &k152, .n1 = 103, .n2 = 103 },
{ .k = &k153, .n1 = 100, .n2 = 6 },
{ .k = &k154, .n1 = 100, .n2 = 100 },
{ .k = &k155, .n1 = 100, .n2 = 100 },
{ .k = &k156, .n1 = 100, .n2 = 100 },
{ .k = &k157, .n1 = 102, .n2 = 102 },
{ .k = &k158, .n1 = 102, .n2 = 105 },
{ .k = &k159, .n1 = 4, .n2 = 105 },
{ .k = &k160, .n1 = 5, .n2 = 102 },
{ .k = &k161, .n1 = 100, .n2 = 105 },
{ .k = &k162, .n1 = 100, .n2 = 6 },
{ .k = &k163, .n1 = 105, .n2 = 103 },
{ .k = &k164, .n1 = 12, .n2 = 103 }
};


void tscaleoutput(FILE *fp, double *xp, double **yp, int kount) {
	for (int i = 1; i <= kount; i++) {
		double myr = xp[i]/(60*60*24*365*1e6);
		fprintf(fp, "%G,", myr);
		
		double *y = yp[i];
		for (int j = 0; j < 164; j++) {
			struct reaction r = reactions[j];
			
			double y1 = loadY(y, r.n1), 
				   y2 = loadY(y, r.n2),
				   kv = r.k(TEMP_INIT, y);
				   
			if (kv != 0) fprintf(fp, "%G,", (kv*y1*y2));
			else fprintf(fp, "0,");
		}
		
		// Special case for 165
		fprintf(fp, "%G", (k165(TEMP_INIT, GRAIN_TEMP_INIT)*getnH(y)*getnH(y)));
		
		// End of step output
		fprintf(fp, "\n");
	}
}

double loadY(double *y, int idx) {
	switch (idx) {
		case 100: return getnH(y);
		case 101: return getnHe(y);
		case 102: return getnC(y);
		case 103: return getnO(y);
		case 104: return getne(y);
		case 105: return getnC(y) + getnO(y) + nSi;
		default:
			if (idx >= 1 && idx <= 27) {
				return y[idx];
			} else {
				printf("Unknown y idx %d\n", idx);
				return 0;
			}
	}
}