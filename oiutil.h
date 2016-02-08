#ifndef _OIUTIL_H_
#define _OIUTIL_H_

double **d,*x;

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nvar);
void lubksb(double **a, int n, int *indx, double b[]);

void ludcmp(double **a, int n, int *indx, double *d);

#endif