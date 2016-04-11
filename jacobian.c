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
	double h2 = h*2;
	for (int i = 1; i <= nderv; i++) {
		double (*deriv)(double[]) = derivs_arr[i-1];
		
		// Apply finite step for each variable
		for (int j = 1; j <= nvar; j++) {
			double original = y[j];
			
			y[j] = original + h;
			double forward = deriv(y);
			y[j] = original - h;
			double backward = deriv(y);
			y[j] = original;
			
			dfdy[i][j] = (forward - backward) / h2;
		}
	}
	
	// free_vector(yh, 1, nvar);
}