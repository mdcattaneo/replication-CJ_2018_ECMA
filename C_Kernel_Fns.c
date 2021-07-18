/*
C FUNCTIONS - Kernel Based Estimators
TO COMPILE: gcc -shared -fPIC -o C_Kernel_Fns.so C_Kernel_Fns.c -lm
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <limits.h>
# include <float.h>

# define EPS DBL_EPSILON

double phi(double x) {return(exp(-0.5*pow(x,2))/sqrt(2.*M_PI));};
double Dphi(double x) {return(-x*phi(x));};

double K(double x, int P){
    double Pol = 0.;
    switch (P) {
    case 2: 
        Pol=1.; break;
    case 4:
        Pol=(3.-pow(x,2))/2.; break;
    case 6:
        Pol=(pow(x,4)-10.*pow(x,2)+15.)/8.; break;
    case 8:
        Pol=(-pow(x,6)+21.*pow(x,4)-105.*pow(x,2)+105.)/48.; break;
    case 10:
        Pol=(pow(x,8)-36.*pow(x,6)+378.*pow(x,4)-1260.*pow(x,2)+945.)/384.; break;
    };
    return(Pol * phi(x));
};

int v_ij(int i, int j, int n) {
    return(i*n-(i*(i+1))/2+j-i-1);
};

void Kest(double *theta_hat,
          double *K_ij,
          double *f,
          double *y, double *x, int *d, int *P, int *n, int *B,
          double *h, int *ids) { 

    int i, j, k, b, bi, bj; double temp=0.;

    double K_ii = pow(K(0.,P[0]),d[0])/pow(h[0],d[0]);

    for (i=0; i<n[0]; i++) for (j=i+1; j<n[0]; j++) {
            temp = 1.; for (k=1; k<d[0]; k++) {temp *= K((x[k*n[0]+i]-x[k*n[0]+j])/h[0],P[0]);}
            K_ij[v_ij(i,j,n[0])] = K((x[i]-x[j])/h[0],P[0])*temp/pow(h[0],d[0]);
    }

    for (b=0; b<(B[0]+1); b++) {
        for (i=0; i<n[0]; i++) {
            f[i]=0.;

            for (j=0; j<n[0]; j++) {
                bi = ids[b*n[0]+i]; bj = ids[b*n[0]+j];

                if (bi <  bj) f[i] += K_ij[v_ij(bi,bj,n[0])];
                if (bi == bj) f[i] += K_ii;
                if (bi >  bj) f[i] += K_ij[v_ij(bj,bi,n[0])];
            }
            theta_hat[b] += (y[bi] >= f[i]/n[0]);
        }
        theta_hat[b] /= n[0];
    }
};

void fi(double *f,
         double *x, int *d, int *P, int *n, double *h) { 

    int i, j, k; double temp=0.;

    for (i=0; i<n[0]; i++) {
        f[i] = 0.;
        for (j=0; j<n[0]; j++) {
            temp = 1.; for (k=1; k<d[0]; k++) {temp *= K((x[k*n[0]+i]-x[k*n[0]+j])/h[0],P[0]);}
            f[i] += K((x[i]-x[j])/h[0],P[0])*temp/pow(h[0],d[0]);
        }
        f[i] /= n[0];
    }
};




