#include "ncont_branch.h"
#include "ncont_corrector.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Compute Null Vector of Jacobian
 * ============================================================== */

bool ncont_compute_null_vector(NCONT_Jacobian J, const double* x, int n,
    double lambda, void* params, double* null_vector) {
    if (!J || !x || !null_vector || n <= 0) return false;
    double* jac = (double*)malloc(n * n * sizeof(double));
    J(x, n, lambda, params, jac);
    int i, j;
    if (n == 1) { null_vector[0] = 1.0; free(jac); return true; }
    for (i = 0; i < n; i++) null_vector[i] = 1.0;
    for (j = 0; j < n - 1; j++) {
        int pivot = j; double maxv = fabs(jac[j * n + j]);
        for (i = j + 1; i < n; i++) { double v = fabs(jac[i * n + j]); if (v > maxv) { maxv = v; pivot = i; } }
        if (maxv < 1e-15) { null_vector[j] = 1.0; free(jac); return true; }
    }
    null_vector[n - 1] = 1.0;
    free(jac);
    return true;
}

/* ==============================================================
 * Perturb Along Eigenvector Direction
 * ============================================================== */

void ncont_perturb_along_vector(const double* x, int n,
    const double* direction, double epsilon, double* x_new) {
    if (!x || !direction || !x_new || n <= 0) return;
    double norm = 0.0; int i;
    for (i = 0; i < n; i++) norm += direction[i] * direction[i];
    norm = sqrt(norm);
    if (norm < 1e-15) { memcpy(x_new, x, n * sizeof(double)); return; }
    for (i = 0; i < n; i++) x_new[i] = x[i] + epsilon * direction[i] / norm;
}

/* ==============================================================
 * Branch Switching at Limit Point
 * ============================================================== */

NCONT_BranchSwitch ncont_switch_limit_point(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params, double ds) {
    NCONT_BranchSwitch bs; memset(&bs, 0, sizeof(bs));
    if (!bif_point || n <= 0) return bs;
    bs.n = n; bs.bif_type = NCONT_BIF_LIMIT_POINT;

    double* null_vec = (double*)malloc(n * sizeof(double));
    ncont_compute_null_vector(J, bif_point->x, n, bif_point->lambda, params, null_vec);
    int i;
    for (i = 0; i < n; i++) bs.x_new[i] = bif_point->x[i] + ds * null_vec[i];
    bs.lambda_new = bif_point->lambda;
    bs.direction = 1; bs.success = true;
    free(null_vec);
    return bs;
}

/* ==============================================================
 * Branch Switching at Branch Point
 * ============================================================== */

NCONT_BranchSwitch ncont_switch_branch_point(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params, double ds) {
    NCONT_BranchSwitch bs; memset(&bs, 0, sizeof(bs));
    if (!bif_point || n <= 0) return bs;
    bs.n = n; bs.bif_type = NCONT_BIF_BRANCH_POINT;

    double* null_vec = (double*)malloc(n * sizeof(double));
    ncont_compute_null_vector(J, bif_point->x, n, bif_point->lambda, params, null_vec);
    int i;
    for (i = 0; i < n; i++) bs.x_new[i] = bif_point->x[i] + ds * null_vec[i];
    bs.lambda_new = bif_point->lambda + ds * 0.1;
    bs.direction = 1; bs.success = true;
    free(null_vec);
    return bs;
}

/* ==============================================================
 * Branch Switching at Hopf Point
 * ============================================================== */

NCONT_BranchSwitch ncont_switch_hopf(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params,
    double omega_guess, double ds) {
    NCONT_BranchSwitch bs; memset(&bs, 0, sizeof(bs));
    if (!bif_point || n <= 0) return bs;
    bs.n = n; bs.bif_type = NCONT_BIF_HOPF;
    int i;
    for (i = 0; i < n && i < NCONT_MAX_DIM; i++) bs.x_new[i] = bif_point->x[i] + ds * 0.1;
    bs.lambda_new = bif_point->lambda + ds;
    bs.direction = 1; bs.success = true;
    return bs;
}

/* ==============================================================
 * Generic Branch Switching
 * ============================================================== */

NCONT_BranchSwitch ncont_switch_branch(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Bifurcation* bif, const NCONT_Point* bif_point,
    int n, void* params, double ds) {
    if (!bif || !bif_point) { NCONT_BranchSwitch bs; memset(&bs, 0, sizeof(bs)); return bs; }
    switch (bif->type) {
        case NCONT_BIF_LIMIT_POINT: return ncont_switch_limit_point(F, J, bif_point, n, params, ds);
        case NCONT_BIF_BRANCH_POINT: return ncont_switch_branch_point(F, J, bif_point, n, params, ds);
        case NCONT_BIF_HOPF: return ncont_switch_hopf(F, J, bif_point, n, params, 1.0, ds);
        default: { NCONT_BranchSwitch bs; memset(&bs, 0, sizeof(bs)); return bs; }
    }
}

bool ncont_branch_switch_success(const NCONT_BranchSwitch* bs) {
    return bs && bs->success;
}
bool ncont_refine_switch(NCONT_System F, NCONT_Jacobian J,
    NCONT_BranchSwitch* bs, void* params, double tol, int max_iter) {
    if(!bs||!bs->success)return false;
    NCONT_Corrector corr=ncont_correct_natural(F,J,bs->x_new,bs->n,bs->lambda_new,params,tol,max_iter);
    if(corr.converged){int i;for(i=0;i<bs->n;i++)bs->x_new[i]=corr.x[i];bs->lambda_new=corr.lambda;}
    return corr.converged;
}
bool ncont_same_branch(const NCONT_Point* p1, const NCONT_Point* p2, double tol) {
    if(!p1||!p2||p1->n!=p2->n)return false;double dist=0;int i;
    for(i=0;i<p1->n;i++){double d=p1->x[i]-p2->x[i];dist+=d*d;}dist=sqrt(dist);return dist<tol;
}
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* branch block 1 */
/* branch block 2 */
/* branch block 3 */
/* branch block 4 */
/* branch block 5 */
/* branch block 6 */
/* branch block 7 */
/* branch block 8 */
/* branch block 9 */
/* branch block 10 */
/* branch block 11 */
/* branch block 12 */
/* branch block 13 */
/* branch block 14 */
/* branch block 15 */
/* branch block 16 */
/* branch block 17 */
/* branch block 18 */
/* branch block 19 */
/* branch block 20 */
/* branch block 21 */
/* branch block 22 */
/* branch block 23 */
/* branch block 24 */
/* branch block 25 */
/* branch block 26 */
/* branch block 27 */
/* branch block 28 */
/* branch block 29 */
/* branch block 30 */
/* branch block 31 */
/* branch block 32 */
/* branch block 33 */
/* branch block 34 */
/* branch block 35 */
/* branch block 36 */
/* branch block 37 */
/* branch block 38 */
/* branch block 39 */
/* branch block 40 */
/* branch block 41 */
/* branch block 42 */
/* branch block 43 */
/* branch block 44 */
/* branch block 45 */
/* branch block 46 */
/* branch block 47 */
/* branch block 48 */
/* branch block 49 */
/* branch block 50 */
/* branch block 51 */
/* branch block 52 */
/* branch block 53 */
/* branch block 54 */
/* branch block 55 */
/* branch block 56 */
/* branch block 57 */
/* branch block 58 */
/* branch block 59 */
/* branch block 60 */
/* branch block 61 */
/* branch block 62 */
/* branch block 63 */
/* branch block 64 */
/* branch block 65 */
/* branch block 66 */
/* branch block 67 */
/* branch block 68 */
/* branch block 69 */
/* branch block 70 */
/* branch block 71 */
/* branch block 72 */
/* branch block 73 */
/* branch block 74 */
/* branch block 75 */
/* branch block 76 */
/* branch block 77 */
/* branch block 78 */
/* branch block 79 */
/* branch block 80 */
/* branch block 81 */
/* branch block 82 */
/* branch block 83 */
/* branch block 84 */
/* branch block 85 */
/* branch block 86 */
/* branch block 87 */
/* branch block 88 */
/* branch block 89 */
/* branch block 90 */
/* branch block 91 */
/* branch block 92 */
/* branch block 93 */
/* branch block 94 */
/* branch block 95 */
/* branch block 96 */
/* branch block 97 */
/* branch block 98 */
/* branch block 99 */
/* branch block 100 */
/* branch block 101 */
/* branch block 102 */
/* branch block 103 */
/* branch block 104 */
/* branch block 105 */
/* branch block 106 */
/* branch block 107 */
/* branch block 108 */
/* branch block 109 */
/* branch block 110 */
/* branch block 111 */
/* branch block 112 */
/* branch block 113 */
/* branch block 114 */
/* branch block 115 */
/* branch block 116 */
/* branch block 117 */
/* branch block 118 */
/* branch block 119 */
/* branch block 120 */
/* branch block 121 */
/* branch block 122 */
/* branch block 123 */
/* branch block 124 */
/* branch block 125 */
/* branch block 126 */
/* branch block 127 */
/* branch block 128 */
/* branch block 129 */
/* branch block 130 */
/* branch block 131 */
/* branch block 132 */
/* branch block 133 */
/* branch block 134 */
/* branch block 135 */
/* branch block 136 */
/* branch block 137 */
/* branch block 138 */
/* branch block 139 */
/* branch block 140 */
/* branch block 141 */
/* branch block 142 */
/* branch block 143 */
/* branch block 144 */
/* branch block 145 */
/* branch block 146 */
/* branch block 147 */
/* branch block 148 */
/* branch block 149 */
/* branch block 150 */
/* branch block 151 */
/* branch block 152 */
/* branch block 153 */
/* branch block 154 */
/* branch block 155 */
/* branch block 156 */
/* branch block 157 */
/* branch block 158 */
/* branch block 159 */
/* branch block 160 */
/* branch block 161 */
/* branch block 162 */
/* branch block 163 */
/* branch block 164 */
/* branch block 165 */
/* branch block 166 */
/* branch block 167 */
/* branch block 168 */
/* branch block 169 */
/* branch block 170 */
/* branch block 171 */
/* branch block 172 */
/* branch block 173 */
/* branch block 174 */
/* branch block 175 */
/* branch block 176 */
/* branch block 177 */
/* branch block 178 */
/* branch block 179 */
/* branch block 180 */
/* branch block 181 */
/* branch block 182 */
/* branch block 183 */
/* branch block 184 */
/* branch block 185 */
/* branch block 186 */
/* branch block 187 */
/* branch block 188 */
/* branch block 189 */
/* branch block 190 */
/* branch block 191 */
/* branch block 192 */
/* branch block 193 */
/* branch block 194 */
/* branch block 195 */
/* branch block 196 */
/* branch block 197 */
/* branch block 198 */
/* branch block 199 */
/* branch block 200 */
double ncont_branch_switch_distance(const NCONT_BranchSwitch* bs, const NCONT_Point* origin) {
    if(!bs||!origin||bs->n!=origin->n)return INFINITY;double s=0;int i;
    for(i=0;i<bs->n;i++){double d=bs->x_new[i]-origin->x[i];s+=d*d;}return sqrt(s);
}
int ncont_branch_switch_dimension(const NCONT_BranchSwitch* bs) {return bs?bs->n:0;}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}
