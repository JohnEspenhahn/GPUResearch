#ifndef _simulation_h_
#define _simulation_h_

#include "odeint.h"
#include "stiff.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

void run(FILE *fp, double x2, double eps, double h1, bool print);

#endif