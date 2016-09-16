// John Espenhahn

#include "jacobian.h"

/**
 * @param derivs Array of derivatives
 * @param nderv Number of derivatives
 * @param y Vector of initial y values
 * @param h Finite step amount
 * @param dfdy Matrix (nderv x nvar) for partials output
 * @param nvar Size of y values vector
 */
void jacobian(double (*derivs_arr[])(double[]), int nderv, double y[], double h, double **dfdy, int nvar) {
	// double *yh = vector(1,nvar);
	// copy_vector(y, yh, 1, nvar);
	
	// Which derivative
	for (int i = 1; i <= nderv; i++) {
		double (*deriv)(double[]) = derivs_arr[i-1];
		double x0 = deriv(y);
		
		// Apply finite step for each variable
		for (int j = 1; j <= nvar; j++) {
			y[j] += h;
			dfdy[i][j] = (deriv(y) - x0) / h;
			y[j] -= h;
			
			/*
			yh[j] += h;
			dfdy[i][j] = (deriv(yh) - x0) / h;
			yh[j] -= h;
			*/
		}
	}
	
	// free_vector(yh, 1, nvar);
}