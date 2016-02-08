#ifndef _STIFBS_H_
#define _STIFBS_H_

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, int, double [], double []));
	
#endif