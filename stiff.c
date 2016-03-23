#include "stiff.h"
#include "oiutil.h"

#define MAXTRY 80 // default 40

#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define TINY 1.0e-30		// Prevents divbyzero
// Here NMAX is the maximum value of n; GROW and SHRNK are the largest and smallest factors
// by which stepsize can change in one step; ERRCON equals (GROW/SAFETY) raised to the power
// (1/PGROW) and handles the case when errmax ~= 0.
#define GAM 0.231
#define A21 2.0
#define A31 4.52470820736
#define A32 4.16352878860
#define C21 -5.07167533877
#define C31 6.02015272865
#define C32 0.159750684673
#define C41 -1.856343618677
#define C42 -8.50538085819
#define C43 -2.08407513602
#define B1 3.95750374663
#define B2 4.62489238836
#define B3 0.617477263873
#define B4 1.282612945268
#define E1 -2.30215540292
#define E2 -3.07363448539
#define E3 0.873280801802
#define E4 1.282612945268

#define C1X GAM
#define C2X -0.396296677520e-01
#define C3X 0.550778939579
#define C4X -0.553509845700e-01
#define A2X 0.462
#define A3X 0.880208333333

void stiff(double y[], double dydx[], int n, double *x, double htry, double eps,
			double yscal[], double *hdid, double *hnext,
			void (*derivs)(double, int, double[], double[]))
// Fourth-order Rosenbrock step for integrating stiffo.d.e.’s, with monitoring of local truncation
// error to adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative
// dydx[1..n] at the starting value of the independent variable x. Also input are the stepsize to
// be attempted htry, the required accuracy eps, and the vector yscal[1..n] against which
// the error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize
// that was actually accomplished, and hnext is the estimated next stepsize. derivs is a usersupplied
// routine that computes the derivatives of the right-hand side with respect to x, while
// jacobn (a fixed name) is a user-supplied routine that computes the Jacobi matrix of derivatives
// of the right-hand side with respect to the components of y.
{
	void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,jtry,*indx;
	double d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err;
	double *g1,*g2,*g3,*g4,*ysav;

	indx=ivector(1,n);
	a=matrix(1,n,1,n);
	dfdx=vector(1,n);
	dfdy=matrix(1,n,1,n);
	dysav=vector(1,n);
	err=vector(1,n);
	g1=vector(1,n);
	g2=vector(1,n);
	g3=vector(1,n);
	g4=vector(1,n);
	ysav=vector(1,n);
	xsav=(*x); // Save initial values.
	for (i=1;i<=n;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	jacobn(xsav,ysav,dfdx,dfdy,n);
	// The user must supply this routine to return the n-by-n matrix dfdy and the vector dfdx.
	h=htry; // Set stepsize to the initial trial value.
	for (jtry=1;jtry<=MAXTRY;jtry++) {
		for (i=1;i<=n;i++) { // Set up the matrix 1 - gamma*h*f'
			for (j=1;j<=n;j++) a[i][j] = -dfdy[i][j];
			a[i][i] += 1.0/(GAM*h);
		}
		ludcmp(a,n,indx,&d); // LU decomposition of the matrix.
		for (i=1;i<=n;i++) // Set up right-hand side for g1.
			g1[i]=dysav[i]+h*C1X*dfdx[i];
		
		lubksb(a,n,indx,g1); // Solve for g1.
		for (i=1;i<=n;i++) // Compute intermediate values of y and x.
			y[i]=ysav[i]+A21*g1[i];
		*x=xsav+A2X*h;
		(*derivs)(*x,n,y,dydx); // Compute dydx at the intermediate values.
		for (i=1;i<=n;i++) // Set up right-hand side for g2.
			g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
		lubksb(a,n,indx,g2); // Solve for g2.
		for (i=1;i<=n;i++) // Compute intermediate values of y and x.
			y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
			
		*x=xsav+A3X*h;
		(*derivs)(*x,n,y,dydx); // Compute dydx at the intermediate values.
		for (i=1;i<=n;i++) // Set up right-hand side for g3.
			g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
		lubksb(a,n,indx,g3); // Solve for g3.
		for (i=1;i<=n;i++) // Set up right-hand side for g4.
			g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
		lubksb(a,n,indx,g4); // Solve for g4.
		for (i=1;i<=n;i++) { // Get fourth-order estimate of y and error estimate.
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
		}
		*x=xsav+h;
		if (*x == xsav) nrerror("stepsize not significant in stiff");
		errmax=0.0; // Evaluate accuracy.
		for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(err[i]/yscal[i]));
		
		errmax /= eps; // Scale relative to required tolerance.
		if (errmax <= 1.0) { // Step succeeded. Compute size of next step and re-turn
			*hdid=h;
			*hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
			free_vector(ysav,1,n);
			free_vector(g4,1,n);
			free_vector(g3,1,n);
			free_vector(g2,1,n);
			free_vector(g1,1,n);
			free_vector(err,1,n);
			free_vector(dysav,1,n);
			free_matrix(dfdy,1,n,1,n);
			free_vector(dfdx,1,n);
			free_matrix(a,1,n,1,n);
			free_ivector(indx,1,n);
			return;
		} else { // Truncation error too large, reduce step size
			*hnext=SAFETY*h*pow(errmax,PSHRNK);
			h=(h >= 0.0 ? FMAX(*hnext,SHRNK*h) : FMIN(*hnext,SHRNK*h));
		}
	} // Go back and re-try step
	nrerror("exceeded MAXTRY in stiff");
}