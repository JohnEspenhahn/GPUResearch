#ifndef _simulation_h_
#define _simulation_h_

#include "odeint.h"
#include "stiff.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define KMAX 100

void run(FILE *fp, double x2, bool print);

#endif