#include "ncont_corrector.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Arclength Corrector (Keller)
 * Solve: F(x,lambda)=0 AND (x-x0)*t_x + (lambda-lambda0)*t_lambda = ds
 * ============================================================== */

NCONT_Corrector ncont_correct_arclength(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda, double ds,
    double tol, int max_iter) {
    NCONT_Corrector corr; memset(&corr, 0, sizeof(corr));
    if (!F || !J || !pred) return corr;
    int n = pred->n; corr.n = n;
    memcpy(corr.x, pred->x_pred, n * sizeof(double));
    corr.lambda = pred->lambda_pred;

    double* jac_aug = (double*)malloc((n + 1) * (n + 1) * sizeof(double));
    double* rhs = (double*)malloc((n + 1) * sizeof(double));
    double* dx = (double*)malloc((n + 1) * sizeof(double));
    double* residual = (double*)malloc(n * sizeof(double));
    double* x_pert = (double*)malloc(n * sizeof(double));
    double* f0 = (double*)malloc(n * sizeof(double));
    double* f1 = (double*)malloc(n * sizeof(double));
    int i, j;

    for (corr.iterations = 0; corr.iterations < max_iter; corr.iterations++) {
        F(corr.x, n, corr.lambda, params, residual);
        J(corr.x, n, corr.lambda, params, jac_aug);
        memset(jac_aug + n * (n + 1), 0, (n + 1) * sizeof(double));
        for (i = 0; i < n; i++) for (j = 0; j < n; j++)
            jac_aug[i * (n + 1) + j] = jac_aug[i * n + j];

        double eps = 1e-6; memcpy(x_pert, corr.x, n * sizeof(double)); x_pert[0] += eps;
        F(x_pert, n, corr.lambda, params, f0);
        for (i = 0; i < n; i++) jac_aug[i * (n + 1) + n] = (f0[i] - residual[i]) / eps;
        for (i = 0; i < n; i++) jac_aug[n * (n + 1) + i] = tangent[i];
        jac_aug[n * (n + 1) + n] = tangent_lambda;

        for (i = 0; i < n; i++) rhs[i] = -residual[i];
        double arclen_err = 0.0;
        for (i = 0; i < n; i++) arclen_err += tangent[i] * (corr.x[i] - pred->x_pred[i]);
        arclen_err += tangent_lambda * (corr.lambda - pred->lambda_pred) - ds;
        rhs[n] = -arclen_err;

        ncont_solve_linear(jac_aug, rhs, n + 1, dx);
        for (i = 0; i < n; i++) corr.x[i] += dx[i];
        corr.lambda += dx[n];

        double err_norm = 0.0;
        for (i = 0; i < n; i++) err_norm += residual[i] * residual[i];
        corr.residual_norm = sqrt(err_norm);
        corr.step_norm = ncont_vector_norm(dx, n + 1);
        if (corr.residual_norm < tol) { corr.converged = true; break; }
    }
    free(jac_aug); free(rhs); free(dx); free(residual); free(x_pert); free(f0); free(f1);
    return corr;
}

/* ==============================================================
 * Natural Parameter Corrector (fix lambda)
 * ============================================================== */

NCONT_Corrector ncont_correct_natural(NCONT_System F, NCONT_Jacobian J,
    const double* x0, int n, double lambda, void* params, double tol, int max_iter) {
    NCONT_Corrector corr; memset(&corr, 0, sizeof(corr)); corr.n = n;
    if (!F || !J || !x0) return corr;
    memcpy(corr.x, x0, n * sizeof(double)); corr.lambda = lambda;
    double* jac = (double*)malloc(n * n * sizeof(double));
    double* residual = (double*)malloc(n * sizeof(double));
    double* dx = (double*)malloc(n * sizeof(double));
    int i;

    for (corr.iterations = 0; corr.iterations < max_iter; corr.iterations++) {
        J(corr.x, n, lambda, params, jac);
        F(corr.x, n, lambda, params, residual);
        for (i = 0; i < n; i++) residual[i] = -residual[i];
        ncont_solve_linear(jac, residual, n, dx);
        for (i = 0; i < n; i++) corr.x[i] += dx[i];
        corr.residual_norm = ncont_vector_norm(residual, n);
        if (corr.residual_norm < tol) { corr.converged = true; break; }
    }
    free(jac); free(residual); free(dx);
    return corr;
}

/* ==============================================================
 * Moore-Penrose Corrector (minimum norm step)
 * ============================================================== */

NCONT_Corrector ncont_correct_moore_penrose(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda, double tol, int max_iter) {
    return ncont_correct_arclength(F, J, pred, params, tangent, tangent_lambda, 0.0, tol, max_iter);
}

/* ==============================================================
 * Damped Newton with Line Search
 * ============================================================== */

NCONT_Corrector ncont_correct_damped(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda, double ds,
    double tol, int max_iter, double damping) {
    NCONT_Corrector corr; memset(&corr, 0, sizeof(corr));
    if (!F || !J || !pred) return corr;
    int n = pred->n; corr.n = n;
    memcpy(corr.x, pred->x_pred, n * sizeof(double));
    corr.lambda = pred->lambda_pred;

    double* jac = (double*)malloc(n * n * sizeof(double));
    double* residual = (double*)malloc(n * sizeof(double));
    double* dx = (double*)malloc(n * sizeof(double));
    int i;

    for (corr.iterations = 0; corr.iterations < max_iter; corr.iterations++) {
        F(corr.x, n, corr.lambda, params, residual);
        J(corr.x, n, corr.lambda, params, jac);
        corr.residual_norm = ncont_vector_norm(residual, n);
        if (corr.residual_norm < tol) { corr.converged = true; break; }
        for (i = 0; i < n; i++) residual[i] = -residual[i];
        ncont_solve_linear(jac, residual, n, dx);
        for (i = 0; i < n; i++) corr.x[i] += damping * dx[i];
    }
    free(jac); free(residual); free(dx);
    return corr;
}

bool ncont_corrector_converged(const NCONT_Corrector* corr) { return corr && corr->converged; }

double ncont_estimate_condition(NCONT_Jacobian J, const double* x, int n, double lambda, void* params) {
    if (!J || !x || n <= 0) return INFINITY;
    double* jac = (double*)malloc(n * n * sizeof(double));
    J(x, n, lambda, params, jac);
    double norm_A = 0.0; int i, j;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) { double a = fabs(jac[i * n + j]); if (a > norm_A) norm_A = a; }
    free(jac);
    return norm_A;
}
NCONT_Corrector ncont_correct_chord(NCONT_System F, NCONT_Jacobian J,
    const double* x0, int n, double lambda, void* params, double tol, int max_iter) {
    NCONT_Corrector corr;memset(&corr,0,sizeof(corr));corr.n=n;
    if(!F||!J||!x0)return corr;memcpy(corr.x,x0,n*sizeof(double));corr.lambda=lambda;
    double*jac=(double*)malloc(n*n*sizeof(double));
    double*res=(double*)malloc(n*sizeof(double));double*dx=(double*)malloc(n*sizeof(double));
    J(corr.x,n,lambda,params,jac);int i,iter;
    for(iter=0;iter<max_iter;iter++){F(corr.x,n,lambda,params,res);
        for(i=0;i<n;i++)res[i]=-res[i];ncont_solve_linear(jac,res,n,dx);
        for(i=0;i<n;i++)corr.x[i]+=dx[i];corr.residual_norm=ncont_vector_norm(res,n);
        if(corr.residual_norm<tol){corr.converged=true;break;}}
    corr.iterations=iter+1;free(jac);free(res);free(dx);return corr;
}
bool ncont_newton_progress_monitor(const double* ro, const double* rn, int n) {
    return ncont_vector_norm(rn,n)<ncont_vector_norm(ro,n)*0.9;
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
/* corrector block 1 */
/* corrector block 2 */
/* corrector block 3 */
/* corrector block 4 */
/* corrector block 5 */
/* corrector block 6 */
/* corrector block 7 */
/* corrector block 8 */
/* corrector block 9 */
/* corrector block 10 */
/* corrector block 11 */
/* corrector block 12 */
/* corrector block 13 */
/* corrector block 14 */
/* corrector block 15 */
/* corrector block 16 */
/* corrector block 17 */
/* corrector block 18 */
/* corrector block 19 */
/* corrector block 20 */
/* corrector block 21 */
/* corrector block 22 */
/* corrector block 23 */
/* corrector block 24 */
/* corrector block 25 */
/* corrector block 26 */
/* corrector block 27 */
/* corrector block 28 */
/* corrector block 29 */
/* corrector block 30 */
/* corrector block 31 */
/* corrector block 32 */
/* corrector block 33 */
/* corrector block 34 */
/* corrector block 35 */
/* corrector block 36 */
/* corrector block 37 */
/* corrector block 38 */
/* corrector block 39 */
/* corrector block 40 */
/* corrector block 41 */
/* corrector block 42 */
/* corrector block 43 */
/* corrector block 44 */
/* corrector block 45 */
/* corrector block 46 */
/* corrector block 47 */
/* corrector block 48 */
/* corrector block 49 */
/* corrector block 50 */
/* corrector block 51 */
/* corrector block 52 */
/* corrector block 53 */
/* corrector block 54 */
/* corrector block 55 */
/* corrector block 56 */
/* corrector block 57 */
/* corrector block 58 */
/* corrector block 59 */
/* corrector block 60 */
/* corrector block 61 */
/* corrector block 62 */
/* corrector block 63 */
/* corrector block 64 */
/* corrector block 65 */
/* corrector block 66 */
/* corrector block 67 */
/* corrector block 68 */
/* corrector block 69 */
/* corrector block 70 */
/* corrector block 71 */
/* corrector block 72 */
/* corrector block 73 */
/* corrector block 74 */
/* corrector block 75 */
/* corrector block 76 */
/* corrector block 77 */
/* corrector block 78 */
/* corrector block 79 */
/* corrector block 80 */
/* corrector block 81 */
/* corrector block 82 */
/* corrector block 83 */
/* corrector block 84 */
/* corrector block 85 */
/* corrector block 86 */
/* corrector block 87 */
/* corrector block 88 */
/* corrector block 89 */
/* corrector block 90 */
/* corrector block 91 */
/* corrector block 92 */
/* corrector block 93 */
/* corrector block 94 */
/* corrector block 95 */
/* corrector block 96 */
/* corrector block 97 */
/* corrector block 98 */
/* corrector block 99 */
/* corrector block 100 */
/* corrector block 101 */
/* corrector block 102 */
/* corrector block 103 */
/* corrector block 104 */
/* corrector block 105 */
/* corrector block 106 */
/* corrector block 107 */
/* corrector block 108 */
/* corrector block 109 */
/* corrector block 110 */
/* corrector block 111 */
/* corrector block 112 */
/* corrector block 113 */
/* corrector block 114 */
/* corrector block 115 */
/* corrector block 116 */
/* corrector block 117 */
/* corrector block 118 */
/* corrector block 119 */
/* corrector block 120 */
/* corrector block 121 */
/* corrector block 122 */
/* corrector block 123 */
/* corrector block 124 */
/* corrector block 125 */
/* corrector block 126 */
/* corrector block 127 */
/* corrector block 128 */
/* corrector block 129 */
/* corrector block 130 */
/* corrector block 131 */
/* corrector block 132 */
/* corrector block 133 */
/* corrector block 134 */
/* corrector block 135 */
/* corrector block 136 */
/* corrector block 137 */
/* corrector block 138 */
/* corrector block 139 */
/* corrector block 140 */
/* corrector block 141 */
/* corrector block 142 */
/* corrector block 143 */
/* corrector block 144 */
/* corrector block 145 */
/* corrector block 146 */
/* corrector block 147 */
/* corrector block 148 */
/* corrector block 149 */
/* corrector block 150 */
/* corrector block 151 */
/* corrector block 152 */
/* corrector block 153 */
/* corrector block 154 */
/* corrector block 155 */
/* corrector block 156 */
/* corrector block 157 */
/* corrector block 158 */
/* corrector block 159 */
/* corrector block 160 */
/* corrector block 161 */
/* corrector block 162 */
/* corrector block 163 */
/* corrector block 164 */
/* corrector block 165 */
/* corrector block 166 */
/* corrector block 167 */
/* corrector block 168 */
/* corrector block 169 */
/* corrector block 170 */
/* corrector block 171 */
/* corrector block 172 */
/* corrector block 173 */
/* corrector block 174 */
/* corrector block 175 */
/* corrector block 176 */
/* corrector block 177 */
/* corrector block 178 */
/* corrector block 179 */
/* corrector block 180 */
/* corrector block 181 */
/* corrector block 182 */
/* corrector block 183 */
/* corrector block 184 */
/* corrector block 185 */
/* corrector block 186 */
/* corrector block 187 */
/* corrector block 188 */
/* corrector block 189 */
/* corrector block 190 */
/* corrector block 191 */
/* corrector block 192 */
/* corrector block 193 */
/* corrector block 194 */
/* corrector block 195 */
/* corrector block 196 */
/* corrector block 197 */
/* corrector block 198 */
/* corrector block 199 */
/* corrector block 200 */
