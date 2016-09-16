// Numerical Recipes in C, The Art of Scientific Computing, Second Edition

#include "oiutil.h"
#include "nrutil.h"
#include <math.h>

#define TINY 1.0e-30		// Prevents divbyzero

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

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}

#undef TINY