#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include "nrutil.h"

void jacobian(double (*derivs[])(double[]), int nderv, double y[], double h, double **dfdy, int nvar);

#endif