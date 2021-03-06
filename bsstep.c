// Numerical Recipes in C, The Art of Scientific Computing, Second Edition

#include "bsstep.h"

#define KMAXX 16			// Max rows in extrapolation. Default = 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25			// Safety factors
#define SAFE2 0.7
#define REDMAX 1.0e-5		// Max factor for stepsize reduction
#define REDMIN 0.7			// Min factor for stepsize reduction
#define SCALMX 0.1			// 1/SCALMX is max factor by which stepsize can increate
#define TINY 1.0e-30		// Prevents divbyzero

double **d, *x; // Pointers to matrix and vector used by pzextr

void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[], void (*derivs)(double, int, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);

void bsstep(double y[], double dydx[], int nvar, double *xx, double htry, 
		double eps, double yscal[], double *hdid, double *hnext,
		void (*derivs)(double, int, double[], double[])) {
// Bulirsch-Stoer step with monitoring of local truncation error to 
// ensure accuracy and adjust stepsize. Input are the dependent 
// variable vector y[1..nvar] and its derivative dydx[1..nvar] at the 
// starting value of the independent variable x. Also input are the 
// stepsize to be attempted htry, the required accuracy eps, and the 
// vector yscal[1..nvar] against which the error is scaled. On output, 
// y and x are replaced by their new values, hdid is the stepsize that
// was actually accomplished, and hnext is the estimated next stepsize.
// derivs is the user-supplied routine that computes the right-hand 
// side derivatives. Be sure to set htry on successive steps to the 
// value of hnext returned from the previous step, as is the case if 
// the routine is called by odeint.
	int i,iq,k,kk,km; 
	static int first=1,kmax,kopt;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *err,*yerr,*ysav,*yseq;
	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
	int reduct,exitflag=0;
	d=matrix(1,nvar,1,KMAXX);
	err=vector(1,KMAXX);
	x=vector(1,KMAXX);
	yerr=vector(1,nvar);
	ysav=vector(1,nvar);
	yseq=vector(1,nvar);
	if (eps != epsold) { // A new tolerance, so reinitialize. 
		*hnext = xnew = -1.0e29; // “Impossible” values. 
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1; // Compute work coefcients Ak. 
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		
		for (iq=2;iq<=KMAXX;iq++) { // Compute a(k,q). 
			for (k=1;k<iq;k++) 
				alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.0)*(2*k+1)));
		} 
		epsold=eps;
		for (kopt=2;kopt<KMAXX;kopt++) // Determine optimal row number for convergence.
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	} 
	h=htry;
	for (i=1;i<=nvar;i++) ysav[i]=y[i]; // Save the starting values. 
	if (*xx != xnew || h != (*hnext)) { // A new stepsize or a new integration: 
										// re-establish the order window.
		first=1;
		kopt=kmax;
	} 
	reduct=0;
	for (;;) { 
		for (k=1;k<=kmax;k++) { // Evaluate the sequence of modi?ed midpoint integrations.
			xnew=(*xx)+h;
			if (xnew == (*xx)) nrerror("step size underflow in bsstep");
			mmid(ysav,dydx,nvar,*xx,h,nseq[k],yseq,derivs);
			xest=SQR(h/nseq[k]); // Squared, since error series is even. 
			pzextr(k,xest,yseq,y,yerr,nvar); // Perform extrapolation. 
			if (k != 1) { // Compute normalized error estimate sigma(k).
				errmax=TINY;
				for (i=1;i<=nvar;i++) 
					errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps; // Scale error relative to tolerance. 
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			} 
			if (k != 1 && (k >= kopt-1 || first)) { // In order window. 
				if (errmax < 1.0) { // Converged. 
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) { // Check for possible stepsize reduction.
					red=SAFE2/err[km];
					break;
				} else if (k == kopt && alf[kopt-1][kopt] < err[km]) { 
					red=1.0/err[km];
					break;
				} else if (kopt == kmax && alf[km][kmax-1] < err[km]) { 
					red=alf[km][kmax-1]*SAFE2/err[km];
					break;
				} else if (alf[km][kopt] < err[km]) { 
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		} 
		if (exitflag) break;
		red=FMIN(red,REDMIN); // Reduce stepsize by at least REDMIN and at most REDMAX.
		red=FMAX(red,REDMAX);
		h *= red;
		reduct=1;
	} // Try again. 
	*xx=xnew; // Successful step taken. 
	*hdid=h;
	first=0;
	wrkmin=1.0e35; // Compute optimal row for convergence and corresponding stepsize.
	for (kk=1;kk<=km;kk++) { 
		fact=FMAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) { 
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		} 
	} 
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) { 
		// Check for possible order increase, but not if stepsize was just reduced. 
		fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) { 
			*hnext=h/fact;
			kopt++;
		} 
	} 
	free_vector(yseq,1,nvar);
	free_vector(ysav,1,nvar);
	free_vector(yerr,1,nvar);
	free_vector(x,1,KMAXX);
	free_vector(err,1,KMAXX);
	free_matrix(d,1,nvar,1,KMAXX);
}

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar) 
// Use polynomial extrapolation to evaluate nvar functions at x = 0 by
// ?tting a polynomial to a sequence of estimates with progressively 
// smaller values x = xest, and corresponding function vectors 
// yest[1..nvar]. This call is number iest in the sequence of calls. 
// Extrapolated function values are output as yz[1..nvar], and their 
// estimated error is output as dy[1..nvar]. 
{ 
	int k1,j; 
	double q,f2,f1,delta,*c;
	
	c=vector(1,nvar); 
	x[iest]=xest; // Save current independent variable. 
	for (j=1;j<=nvar;j++) dy[j]=yz[j]=yest[j]; 
	
	if (iest == 1) { // Store ?rst estimate in ?rst column. 
		for (j=1;j<=nvar;j++) d[j][1]=yest[j]; 
	} else { 
		for (j=1;j<=nvar;j++) c[j]=yest[j]; 
		for (k1=1;k1<iest;k1++) { 
			delta=1.0/(x[iest-k1]-xest); 
			f1=xest*delta; 
			f2=x[iest-k1]*delta; 
			for (j=1;j<=nvar;j++) { // Propagate tableau 1 diagonal more. 
				q=d[j][k1]; 
				d[j][k1]=dy[j]; 
				delta=c[j]-q;
				dy[j]=f1*delta; 
				c[j]=f2*delta; 
				yz[j] += dy[j]; 
			} 
		} 
		for (j=1;j<=nvar;j++) 
			d[j][iest]=dy[j]; 
	} 
	free_vector(c,1,nvar);
}


void mmid(double y[], double dydx[], int nvar, double xs, double htot, 
		int nstep, double yout[], void (*derivs)(double, int, double[], double[])) 
// Modifed midpoint step. At xs, input the dependent variable vector y[1..nvar] and its derivative vector dydx[1..nvar]. Also input is htot, 
// the total step to be made, and nstep, the number of substeps to be used. The output is returned as yout[1..nvar], which need not be a distinct
// array from y; if it is distinct, however, then y and dydx are returned undamaged. 
{ 
	int n,i; 
	double x,swap,h2,h,*ym,*yn;

	ym=vector(1,nvar); 
	yn=vector(1,nvar); 
	h=htot/nstep; // Stepsize this trip
	
	for (i=1;i<=nvar;i++) { 
		ym[i]=y[i]; 
		yn[i]=y[i]+h*dydx[i]; // First step. 
	}

	x = xs + h;
	(*derivs)(x,nvar,yn,yout); // Use yout for temp storage of derivatives
	h2=2.0*h;
	for (n = 2; n <= nstep; n++) {
		for (i = 1; i <= nvar; i++) {
			swap = ym[i]+h2*yout[i];
			ym[i] = yn[i];
			yn[i] = swap;
		}
		x += h;
		(*derivs)(x,nvar,yn,yout);
	}

	for (i = 1; i <= nvar; i++) // Last step
		yout[i] = 0.5 * (ym[i]+yn[i]+h*yout[i]);

	free_vector(yn,1,nvar);
	free_vector(ym,1,nvar);
}

#undef TINY