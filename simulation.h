#ifndef _simulation_h_
#define _simulation_h_

#include "odeint.h"
#include "derivs.h"
#include "stiff.h"
#include <stdlib.h>
#include <stdio.h>

#define KMAX 100

void run(FILE *fp, double x2);

#endif