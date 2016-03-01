#ifndef _simulation_h_
#define _simulation_h_

#include "odeint.h"
#include "gml/derivs.h"
#include "stiff.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define ZONES 1
#define NVAR 2*ZONES

#define KMAX 100

void run(FILE *fp, double x2, bool print);

#endif