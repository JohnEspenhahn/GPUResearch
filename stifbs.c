// Numerical Recipes in C, The Art of Scientific Computing, Second Edition

#include "nrutil.h"
#include <math.h>
#include "stiff.h"
#include "oiutil.h"

#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, int, double [], double []))
// Semi-implicit extrapolation step for integrating stiffo.d.e.’s, with monitoring of local truncation
// error to adjust stepsize. Input are the dependent variable vector y[1..nv] and its derivative
// dydx[1..nv] at the starting value of the independent variable x. Also input are the stepsize
// to be attempted htry, the required accuracy eps, and the vector yscal[1..nv] against
// which the error is scaled. On output, y and x are replaced by their new values, hdid is the
// stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
// is a user-supplied routine that computes the derivatives of the right-hand side with respect to
// x, while jacobn (a fixed name) is a user-supplied routine that computes the Jacobi matrix of
// derivatives of the right-hand side with respect to the components of y. Be sure to set htry
// on successive steps to the value of hnext returned from the previous step, as is the case if the
// routine is called by odeint.
{
	void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
	void simpr(double y[], double dydx[], double dfdx[], double **dfdy,
		int n, double xs, double htot, int nstep, double yout[],
		void (*derivs)(double, int, double [], double []));
	void pzextr(int iest, double xest, double yest[], double yz[], double dy[],int nv);
	
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70}; // Sequence is different from bsstep.
	int reduct,exitflag=0;
	d=matrix(1,nv,1,KMAXX);
	dfdx=vector(1,nv);
	dfdy=matrix(1,nv,1,nv);
	err=vector(1,KMAXX);
	x=vector(1,KMAXX);
	yerr=vector(1,nv);
	ysav=vector(1,nv);
	yseq=vector(1,nv);
	if(eps != epsold || nv != nvold) { // Reinitialize also if nv has changed.
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold=eps;
		nvold=nv; // Save nv.
		a[1] += nv; // Add cost of Jacobian evaluations to work
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1]; // coefficients.
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) ysav[i]=y[i];
	jacobn(*xx,y,dfdx,dfdy,nv); // Evaluate Jacobian.
	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			if (xnew == (*xx)) nrerror("step size underflow in stifbs");
			simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
			// Semi-implicit midpoint rule.
			xest=SQR(h/nseq[k]);
			pzextr(k,xest,yseq,y,yerr,nv);
			if (k != 1) {
				errmax=TINY;
				for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
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
		red=FMIN(red,REDMIN);
		red=FMAX(red,REDMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
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
		fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
	free_vector(yseq,1,nv);
	free_vector(ysav,1,nv);
	free_vector(yerr,1,nv);
	free_vector(x,1,KMAXX);
	free_vector(err,1,KMAXX);
	free_matrix(dfdy,1,nv,1,nv);
	free_vector(dfdx,1,nv);
	free_matrix(d,1,nv,1,KMAXX);
}

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, int, double [], double []))
// Performs one step of semi-implicit midpoint rule. Input are the dependent variable y[1..n], its
// derivative dydx[1..n], the derivative of the right-hand side with respect to x, dfdx[1..n],
// and the Jacobian dfdy[1..n][1..n] at xs. Also input are htot, the total step to be taken,
// and nstep, the number of substeps to be used. The output is returned as yout[1..n].
// derivs is the user-supplied routine that calculates dydx.
{
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;
	indx=ivector(1,n);
	a=matrix(1,n,1,n);
	del=vector(1,n);
	ytemp=vector(1,n);
	h=htot/nstep; // Stepsize this trip.
	for (i=1;i<=n;i++) { // Set up the matrix 1 - hf.
		for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
		++a[i][i];
	}
	ludcmp(a,n,indx,&d); // LU decomposition of the matrix.
	for (i=1;i<=n;i++) // Set up right-hand side for first step. Use yout
		yout[i]=h*(dydx[i]+h*dfdx[i]); // for temporary storage.
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++) // First step.
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	(*derivs)(x,n,ytemp,yout); // Use yout for temporary storage of derivatives.
	(*derivs)(x,n,ytemp,yout); // Use yout for temporary storage of derivatives.
	for (nn=2;nn<=nstep;nn++) { // General step.
		for (i=1;i<=n;i++) // Set up right-hand side for general step.
			yout[i]=h*yout[i]-del[i];
		lubksb(a,n,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		(*derivs)(x,n,ytemp,yout);
	}
	for (i=1;i<=n;i++) // Set up right-hand side for last step.
		yout[i]=h*yout[i]-del[i];
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++) // Take last step.
		yout[i] += ytemp[i];
		
	free_vector(ytemp,1,n);
	free_vector(del,1,n);
	free_matrix(a,1,n,1,n);
	free_ivector(indx,1,n);
}