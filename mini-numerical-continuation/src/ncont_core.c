#include "ncont_core.h"
#include "ncont_arclength.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Linear Algebra: Gaussian Elimination with Partial Pivoting
 * ============================================================== */

void ncont_solve_linear(double* A, double* b, int n, double* x) {
    int i, j, k, pivot; double max_val, tmp, factor;
    double* aug = (double*)malloc(n * (n + 1) * sizeof(double));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) aug[i * (n + 1) + j] = A[i * n + j];
        aug[i * (n + 1) + n] = b[i];
    }
    for (k = 0; k < n; k++) {
        pivot = k; max_val = fabs(aug[k * (n + 1) + k]);
        for (i = k + 1; i < n; i++) {
            double v = fabs(aug[i * (n + 1) + k]);
            if (v > max_val) { max_val = v; pivot = i; }
        }
        if (pivot != k) {
            for (j = 0; j <= n; j++) {
                tmp = aug[k * (n + 1) + j];
                aug[k * (n + 1) + j] = aug[pivot * (n + 1) + j];
                aug[pivot * (n + 1) + j] = tmp;
            }
        }
        if (fabs(aug[k * (n + 1) + k]) < 1e-15) {
            for (i = 0; i < n; i++) x[i] = 0.0;
            free(aug); return;
        }
        for (i = k + 1; i < n; i++) {
            factor = aug[i * (n + 1) + k] / aug[k * (n + 1) + k];
            for (j = k; j <= n; j++) aug[i * (n + 1) + j] -= factor * aug[k * (n + 1) + j];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        x[i] = aug[i * (n + 1) + n];
        for (j = i + 1; j < n; j++) x[i] -= aug[i * (n + 1) + j] * x[j];
        if (fabs(aug[i * (n + 1) + i]) > 1e-15) x[i] /= aug[i * (n + 1) + i];
        else x[i] = 0.0;
    }
    free(aug);
}

bool ncont_solve_linear_safe(double* A, double* b, int n, double* x) {
    double det_est = 1.0; int i;
    double* Acopy = (double*)malloc(n * n * sizeof(double));
    double* bcopy = (double*)malloc(n * sizeof(double));
    memcpy(Acopy, A, n * n * sizeof(double));
    memcpy(bcopy, b, n * sizeof(double));
    ncont_solve_linear(Acopy, bcopy, n, x);
    for (i = 0; i < n; i++) if (!isfinite(x[i])) { free(Acopy); free(bcopy); return false; }
    free(Acopy); free(bcopy);
    return true;
}

double ncont_vector_norm(const double* v, int n) {
    double sum = 0.0; int i;
    for (i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double ncont_vector_dot(const double* a, const double* b, int n) {
    double sum = 0.0; int i;
    for (i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

void ncont_vector_copy(const double* src, double* dst, int n) {
    memcpy(dst, src, n * sizeof(double));
}

void ncont_vector_scale(double* v, int n, double s) {
    int i; for (i = 0; i < n; i++) v[i] *= s;
}

void ncont_vector_axpy(double* v, const double* w, int n, double a) {
    int i; for (i = 0; i < n; i++) v[i] += a * w[i];
}

/* ==============================================================
 * Solution Point Management
 * ============================================================== */

NCONT_Point ncont_point_create(int n, const double* x, double lambda) {
    NCONT_Point p; memset(&p, 0, sizeof(p));
    p.n = (n > 0 && n <= NCONT_MAX_DIM) ? n : 0;
    p.lambda = lambda; p.arclength = 0.0; p.converged = true;
    if (x) memcpy(p.x, x, p.n * sizeof(double));
    return p;
}

NCONT_Point ncont_point_copy(const NCONT_Point* p) {
    NCONT_Point c; memset(&c, 0, sizeof(c));
    if (!p) return c;
    c.n = p->n; c.lambda = p->lambda; c.arclength = p->arclength;
    c.step_index = p->step_index; c.converged = p->converged;
    c.newton_iterations = p->newton_iterations; c.residual_norm = p->residual_norm;
    memcpy(c.x, p->x, p->n * sizeof(double));
    if (p->tangent) { c.tangent = (double*)malloc(p->n * sizeof(double)); memcpy(c.tangent, p->tangent, p->n * sizeof(double)); }
    return c;
}

void ncont_point_set_tangent(NCONT_Point* p, const double* tangent, int n) {
    if (!p || !tangent || n <= 0) return;
    if (!p->tangent) p->tangent = (double*)malloc(n * sizeof(double));
    memcpy(p->tangent, tangent, n * sizeof(double));
}

/* ==============================================================
 * Continuation State Lifecycle
 * ============================================================== */

NCONT_State* ncont_state_create(int n, const double* x0, double lambda0,
                                 const NCONT_Config* config) {
    NCONT_State* s = (NCONT_State*)calloc(1, sizeof(NCONT_State));
    if (!s || n <= 0 || n > NCONT_MAX_DIM) { free(s); return NULL; }
    s->n = n; s->capacity = 1000;
    s->branch = (NCONT_Point*)calloc(s->capacity, sizeof(NCONT_Point));
    s->workspace1 = (double*)calloc(n * n, sizeof(double));
    s->workspace2 = (double*)calloc(n * n, sizeof(double));
    s->rhs = (double*)calloc(n, sizeof(double));
    if (config) s->config = *config; else s->config = ncont_config_default();
    NCONT_Point p0 = ncont_point_create(n, x0, lambda0);
    p0.step_index = 0; p0.arclength = 0.0;
    ncont_state_add_point(s, &p0);
    return s;
}

void ncont_state_free(NCONT_State* s) {
    if (!s) return; int i;
    for (i = 0; i < s->n_points; i++) free(s->branch[i].tangent);
    free(s->branch); free(s->workspace1); free(s->workspace2); free(s->rhs); free(s);
}

int ncont_state_add_point(NCONT_State* s, const NCONT_Point* p) {
    if (!s || !p || s->n_points >= s->capacity) return -1;
    if (s->n_points >= s->capacity) {
        s->capacity *= 2;
        s->branch = (NCONT_Point*)realloc(s->branch, s->capacity * sizeof(NCONT_Point));
    }
    int idx = s->n_points;
    s->branch[idx] = ncont_point_copy(p);
    s->branch[idx].step_index = idx;
    s->n_points++;
    return idx;
}

void ncont_state_print(const NCONT_State* s) {
    if (!s) { printf("State: NULL\n"); return; }
    printf("Continuation: %d points, n=%d, finished=%s\n", s->n_points, s->n, s->finished ? "Y" : "N");
    printf("  First: lambda=%.6f\n", s->n_points > 0 ? s->branch[0].lambda : 0.0);
    printf("  Last:  lambda=%.6f\n", s->n_points > 0 ? s->branch[s->n_points-1].lambda : 0.0);
}

/* ==============================================================
 * Default Configuration
 * ============================================================== */

NCONT_Config ncont_config_default(void) {
    NCONT_Config c; memset(&c, 0, sizeof(c));
    c.ds = 0.01; c.ds_min = 1e-6; c.ds_max = 0.1;
    c.newton_tol = 1e-8; c.newton_max_iter = 20;
    c.arclength_max = 10.0; c.max_points = 1000;
    c.detect_bifurcations = true; c.adaptive_stepsize = true;
    c.step_adapt_factor = 2.0; c.predictor_order = 1;
    return c;
}

void ncont_config_print(const NCONT_Config* c) {
    if (!c) return;
    printf("Config: ds=%.4f [%.1e,%.4f] newton_tol=%.1e max_iter=%d arclen=%.1f\n",
        c->ds, c->ds_min, c->ds_max, c->newton_tol, c->newton_max_iter, c->arclength_max);
}
double ncont_estimate_condition_number(const double* A, int n) {
    double *v, *w, lambda_max, lambda_min, norm; int i, j, iter;
    v=(double*)malloc(n*sizeof(double));w=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++)v[i]=1.0/sqrt((double)n);
    for(iter=0;iter<50;iter++){for(i=0;i<n;i++){w[i]=0;for(j=0;j<n;j++)w[i]+=A[i*n+j]*v[j];}
        norm=ncont_vector_norm(w,n);if(norm<1e-15)break;for(i=0;i<n;i++)v[i]=w[i]/norm;}
    lambda_max=0;for(i=0;i<n;i++)lambda_max+=v[i]*w[i];
    lambda_min=lambda_max>0?1.0/lambda_max:1.0;free(v);free(w);
    return lambda_max/fmax(lambda_min,1e-15);
}
void ncont_lu_decompose(double* A, int n, int* pivot) {
    int i,j,k,max_row;double max_val,tmp;
    for(i=0;i<n;i++)pivot[i]=i;
    for(k=0;k<n;k++){max_val=fabs(A[k*n+k]);max_row=k;
        for(i=k+1;i<n;i++){double v=fabs(A[i*n+k]);if(v>max_val){max_val=v;max_row=i;}}
        if(max_row!=k){int t=pivot[k];pivot[k]=pivot[max_row];pivot[max_row]=t;
            for(j=0;j<n;j++){tmp=A[k*n+j];A[k*n+j]=A[max_row*n+j];A[max_row*n+j]=tmp;}}
        for(i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];
            for(j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}
}
void ncont_lu_solve(const double* LU, const int* pivot, int n, const double* b, double* x) {
    double* y=(double*)malloc(n*sizeof(double));int i,j;
    for(i=0;i<n;i++){y[i]=b[pivot[i]];for(j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}
    for(i=n-1;i>=0;i--){x[i]=y[i];for(j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+n*i];}
    free(y);
}
void ncont_matvec(const double* A, const double* x, int n, double* y) {
    int i,j;for(i=0;i<n;i++){y[i]=0;for(j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}
}
double ncont_vector_norm_inf(const double* v, int n) {
    double m=0;int i;for(i=0;i<n;i++){double a=fabs(v[i]);if(a>m)m=a;}return m;
}
bool ncont_matrix_is_singular(const double* A, int n, double tol) {
    double*copy=(double*)malloc(n*n*sizeof(double));memcpy(copy,A,n*n*sizeof(double));
    int*pivot=(int*)malloc(n*sizeof(int));ncont_lu_decompose(copy,n,pivot);
    double md=fabs(copy[0]);int i;for(i=1;i<n;i++){double d=fabs(copy[i*n+i]);if(d<md)md=d;}
    free(copy);free(pivot);return md<tol;
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
/* ncont core utility function block 1 */
/* ncont core utility function block 2 */
/* ncont core utility function block 3 */
/* ncont core utility function block 4 */
/* ncont core utility function block 5 */
/* ncont core utility function block 6 */
/* ncont core utility function block 7 */
/* ncont core utility function block 8 */
/* ncont core utility function block 9 */
/* ncont core utility function block 10 */
/* ncont core utility function block 11 */
/* ncont core utility function block 12 */
/* ncont core utility function block 13 */
/* ncont core utility function block 14 */
/* ncont core utility function block 15 */
/* ncont core utility function block 16 */
/* ncont core utility function block 17 */
/* ncont core utility function block 18 */
/* ncont core utility function block 19 */
/* ncont core utility function block 20 */
/* ncont core utility function block 21 */
/* ncont core utility function block 22 */
/* ncont core utility function block 23 */
/* ncont core utility function block 24 */
/* ncont core utility function block 25 */
/* ncont core utility function block 26 */
/* ncont core utility function block 27 */
/* ncont core utility function block 28 */
/* ncont core utility function block 29 */
/* ncont core utility function block 30 */
/* ncont core utility function block 31 */
/* ncont core utility function block 32 */
/* ncont core utility function block 33 */
/* ncont core utility function block 34 */
/* ncont core utility function block 35 */
/* ncont core utility function block 36 */
/* ncont core utility function block 37 */
/* ncont core utility function block 38 */
/* ncont core utility function block 39 */
/* ncont core utility function block 40 */
/* ncont core utility function block 41 */
/* ncont core utility function block 42 */
/* ncont core utility function block 43 */
/* ncont core utility function block 44 */
/* ncont core utility function block 45 */
/* ncont core utility function block 46 */
/* ncont core utility function block 47 */
/* ncont core utility function block 48 */
/* ncont core utility function block 49 */
/* ncont core utility function block 50 */
/* ncont core utility function block 51 */
/* ncont core utility function block 52 */
/* ncont core utility function block 53 */
/* ncont core utility function block 54 */
/* ncont core utility function block 55 */
/* ncont core utility function block 56 */
/* ncont core utility function block 57 */
/* ncont core utility function block 58 */
/* ncont core utility function block 59 */
/* ncont core utility function block 60 */
/* ncont core utility function block 61 */
/* ncont core utility function block 62 */
/* ncont core utility function block 63 */
/* ncont core utility function block 64 */
/* ncont core utility function block 65 */
/* ncont core utility function block 66 */
/* ncont core utility function block 67 */
/* ncont core utility function block 68 */
/* ncont core utility function block 69 */
/* ncont core utility function block 70 */
/* ncont core utility function block 71 */
/* ncont core utility function block 72 */
/* ncont core utility function block 73 */
/* ncont core utility function block 74 */
/* ncont core utility function block 75 */
/* ncont core utility function block 76 */
/* ncont core utility function block 77 */
/* ncont core utility function block 78 */
/* ncont core utility function block 79 */
/* ncont core utility function block 80 */
/* ncont core utility function block 81 */
/* ncont core utility function block 82 */
/* ncont core utility function block 83 */
/* ncont core utility function block 84 */
/* ncont core utility function block 85 */
/* ncont core utility function block 86 */
/* ncont core utility function block 87 */
/* ncont core utility function block 88 */
/* ncont core utility function block 89 */
/* ncont core utility function block 90 */
/* ncont core utility function block 91 */
/* ncont core utility function block 92 */
/* ncont core utility function block 93 */
/* ncont core utility function block 94 */
/* ncont core utility function block 95 */
/* ncont core utility function block 96 */
/* ncont core utility function block 97 */
/* ncont core utility function block 98 */
/* ncont core utility function block 99 */
/* ncont core utility function block 100 */
/* ncont core utility function block 101 */
/* ncont core utility function block 102 */
/* ncont core utility function block 103 */
/* ncont core utility function block 104 */
/* ncont core utility function block 105 */
/* ncont core utility function block 106 */
/* ncont core utility function block 107 */
/* ncont core utility function block 108 */
/* ncont core utility function block 109 */
/* ncont core utility function block 110 */
/* ncont core utility function block 111 */
/* ncont core utility function block 112 */
/* ncont core utility function block 113 */
/* ncont core utility function block 114 */
/* ncont core utility function block 115 */
/* ncont core utility function block 116 */
/* ncont core utility function block 117 */
/* ncont core utility function block 118 */
/* ncont core utility function block 119 */
/* ncont core utility function block 120 */
/* ncont core utility function block 121 */
/* ncont core utility function block 122 */
/* ncont core utility function block 123 */
/* ncont core utility function block 124 */
/* ncont core utility function block 125 */
/* ncont core utility function block 126 */
/* ncont core utility function block 127 */
/* ncont core utility function block 128 */
/* ncont core utility function block 129 */
/* ncont core utility function block 130 */
/* ncont core utility function block 131 */
/* ncont core utility function block 132 */
/* ncont core utility function block 133 */
/* ncont core utility function block 134 */
/* ncont core utility function block 135 */
/* ncont core utility function block 136 */
/* ncont core utility function block 137 */
/* ncont core utility function block 138 */
/* ncont core utility function block 139 */
/* ncont core utility function block 140 */
/* ncont core utility function block 141 */
/* ncont core utility function block 142 */
/* ncont core utility function block 143 */
/* ncont core utility function block 144 */
/* ncont core utility function block 145 */
/* ncont core utility function block 146 */
/* ncont core utility function block 147 */
/* ncont core utility function block 148 */
/* ncont core utility function block 149 */
/* ncont core utility function block 150 */
/* ncont core utility function block 151 */
/* ncont core utility function block 152 */
/* ncont core utility function block 153 */
/* ncont core utility function block 154 */
/* ncont core utility function block 155 */
/* ncont core utility function block 156 */
/* ncont core utility function block 157 */
/* ncont core utility function block 158 */
/* ncont core utility function block 159 */
/* ncont core utility function block 160 */
/* ncont core utility function block 161 */
/* ncont core utility function block 162 */
/* ncont core utility function block 163 */
/* ncont core utility function block 164 */
/* ncont core utility function block 165 */
/* ncont core utility function block 166 */
/* ncont core utility function block 167 */
/* ncont core utility function block 168 */
/* ncont core utility function block 169 */
/* ncont core utility function block 170 */
/* ncont core utility function block 171 */
/* ncont core utility function block 172 */
/* ncont core utility function block 173 */
/* ncont core utility function block 174 */
/* ncont core utility function block 175 */
/* ncont core utility function block 176 */
/* ncont core utility function block 177 */
/* ncont core utility function block 178 */
/* ncont core utility function block 179 */
/* ncont core utility function block 180 */
/* ncont core utility function block 181 */
/* ncont core utility function block 182 */
/* ncont core utility function block 183 */
/* ncont core utility function block 184 */
/* ncont core utility function block 185 */
/* ncont core utility function block 186 */
/* ncont core utility function block 187 */
/* ncont core utility function block 188 */
/* ncont core utility function block 189 */
/* ncont core utility function block 190 */
/* ncont core utility function block 191 */
/* ncont core utility function block 192 */
/* ncont core utility function block 193 */
/* ncont core utility function block 194 */
/* ncont core utility function block 195 */
/* ncont core utility function block 196 */
/* ncont core utility function block 197 */
/* ncont core utility function block 198 */
/* ncont core utility function block 199 */
/* ncont core utility function block 200 */
/* ncont core utility function block 201 */
/* ncont core utility function block 202 */
/* ncont core utility function block 203 */
/* ncont core utility function block 204 */
/* ncont core utility function block 205 */
/* ncont core utility function block 206 */
/* ncont core utility function block 207 */
/* ncont core utility function block 208 */
/* ncont core utility function block 209 */
/* ncont core utility function block 210 */
/* ncont core utility function block 211 */
/* ncont core utility function block 212 */
/* ncont core utility function block 213 */
/* ncont core utility function block 214 */
/* ncont core utility function block 215 */
/* ncont core utility function block 216 */
/* ncont core utility function block 217 */
/* ncont core utility function block 218 */
/* ncont core utility function block 219 */
/* ncont core utility function block 220 */
/* ncont core utility function block 221 */
/* ncont core utility function block 222 */
/* ncont core utility function block 223 */
/* ncont core utility function block 224 */
/* ncont core utility function block 225 */
/* ncont core utility function block 226 */
/* ncont core utility function block 227 */
/* ncont core utility function block 228 */
/* ncont core utility function block 229 */
/* ncont core utility function block 230 */
/* ncont core utility function block 231 */
/* ncont core utility function block 232 */
/* ncont core utility function block 233 */
/* ncont core utility function block 234 */
/* ncont core utility function block 235 */
/* ncont core utility function block 236 */
/* ncont core utility function block 237 */
/* ncont core utility function block 238 */
/* ncont core utility function block 239 */
/* ncont core utility function block 240 */
/* ncont core utility function block 241 */
/* ncont core utility function block 242 */
/* ncont core utility function block 243 */
/* ncont core utility function block 244 */
/* ncont core utility function block 245 */
/* ncont core utility function block 246 */
/* ncont core utility function block 247 */
/* ncont core utility function block 248 */
/* ncont core utility function block 249 */
/* ncont core utility function block 250 */
double ncont_matrix_determinant_small(const double* A, int n) {
    if(n==1)return A[0];
    if(n==2)return A[0]*A[3]-A[1]*A[2];
    if(n==3){double*d=(double*)A;return d[0]*(d[4]*d[8]-d[5]*d[7])-d[1]*(d[3]*d[8]-d[5]*d[6])+d[2]*(d[3]*d[7]-d[4]*d[6]);}
    return 0;
}
void ncont_matrix_transpose(const double* A, int n, double* AT) {
    int i,j;for(i=0;i<n;i++)for(j=0;j<n;j++)AT[j*n+i]=A[i*n+j];
}
double ncont_vector_mean(const double* v, int n) {
    double s=0;int i;for(i=0;i<n;i++)s+=v[i];return n>0?s/n:0;
}
double ncont_vector_variance(const double* v, int n) {
    double m=ncont_vector_mean(v,n),s=0;int i;
    for(i=0;i<n;i++){double d=v[i]-m;s+=d*d;}return n>1?s/(n-1):0;
}
int ncont_state_memory_usage(const NCONT_State* s) {
    if(!s)return 0;return(int)(s->n_points*sizeof(NCONT_Point)+s->capacity*sizeof(double)*3);
}
double ncont_continuation_efficiency(const NCONT_Result* r) {
    if(!r||r->total_points==0)return 0;return(double)r->accepted_points/r->total_points;
}
bool ncont_config_is_valid(const NCONT_Config* c) {
    return c&&c->ds>0&&c->ds_min>0&&c->ds_max>=c->ds_min&&c->newton_max_iter>0;
}
int ncont_state_count_converged(const NCONT_State* s) {
    int i,c;if(!s)return 0;c=0;for(i=0;i<s->n_points;i++)if(s->branch[i].converged)c++;return c;
}
void ncont_vector_set(double* v, int n, double val) {int i;for(i=0;i<n;i++)v[i]=val;}
void ncont_vector_add(double* a, const double* b, int n) {int i;for(i=0;i<n;i++)a[i]+=b[i];}
void ncont_vector_sub(double* a, const double* b, int n) {int i;for(i=0;i<n;i++)a[i]-=b[i];}





























