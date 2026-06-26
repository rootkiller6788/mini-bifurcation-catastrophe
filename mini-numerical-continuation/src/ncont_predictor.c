#include "ncont_predictor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Tangent Computation via Bordered Matrix
 * [ J_x   J_lambda ] [ t_x      ]   [ 0 ]
 * [ t_x^T t_lambda  ] [ t_lambda ] = [ 1 ]
 * ============================================================== */

bool ncont_compute_tangent_bordered(NCONT_System F, NCONT_Jacobian J,
    const double* x, int n, double lambda, void* params,
    double* tangent, double* tangent_lambda) {
    if (!F || !J || !x || !tangent || !tangent_lambda || n <= 0) return false;
    double* jac = (double*)malloc(n * n * sizeof(double));
    double* rhs = (double*)malloc((n + 1) * sizeof(double));
    double* aug = (double*)malloc((n + 1) * (n + 1) * sizeof(double));
    double* result = (double*)malloc((n + 1) * sizeof(double));

    J(x, n, lambda, params, jac);
    memset(aug, 0, (n + 1) * (n + 1) * sizeof(double));
    int i, j;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) aug[i * (n + 1) + j] = jac[i * n + j];
    double* dxdl = (double*)malloc(n * sizeof(double));
    double* x_pert = (double*)malloc(n * sizeof(double));
    memcpy(x_pert, x, n * sizeof(double));
    double eps = 1e-6; x_pert[0] += eps;
    double* f0 = (double*)malloc(n * sizeof(double));
    double* f1 = (double*)malloc(n * sizeof(double));
    F(x, n, lambda, params, f0);
    F(x_pert, n, lambda, params, f1);
    for (i = 0; i < n; i++) { aug[i * (n + 1) + n] = (f1[i] - f0[i]) / eps; aug[n * (n + 1) + i] = tangent[i]; }
    for (i = 0; i < n + 1; i++) rhs[i] = 0.0;
    rhs[n] = 1.0;
    ncont_solve_linear(aug, rhs, n + 1, result);
    double norm = 0.0;
    for (i = 0; i < n; i++) norm += result[i] * result[i];
    norm += result[n] * result[n]; norm = sqrt(norm);
    if (norm > 1e-15) { for (i = 0; i < n; i++) tangent[i] = result[i] / norm; *tangent_lambda = result[n] / norm; }
    free(jac); free(rhs); free(aug); free(result); free(dxdl); free(x_pert); free(f0); free(f1);
    return true;
}

void ncont_compute_tangent(NCONT_System F, NCONT_Jacobian J,
    const double* x, int n, double lambda, void* params,
    double* tangent, double* tangent_lambda) {
    ncont_compute_tangent_bordered(F, J, x, n, lambda, params, tangent, tangent_lambda);
}

/* ==============================================================
 * Tangent Predictor
 * ============================================================== */

NCONT_Predictor ncont_predict_tangent(const NCONT_State* state, int idx, double ds) {
    NCONT_Predictor pred; memset(&pred, 0, sizeof(pred));
    if (!state || idx < 0 || idx >= state->n_points) return pred;
    int n = state->n;
    NCONT_Point* p = &state->branch[idx];
    pred.n = n;
    if (p->tangent) {
        double norm = 0.0; int i;
        for (i = 0; i < n; i++) norm += p->tangent[i] * p->tangent[i];
        if (norm > 1e-15) {
            double scale = ds / sqrt(norm);
            for (i = 0; i < n; i++) pred.x_pred[i] = p->x[i] + scale * p->tangent[i];
            pred.lambda_pred = p->lambda + ds * 0.1;
        }
    }
    pred.predictor_norm = ncont_vector_norm(pred.x_pred, n);
    return pred;
}

/* ==============================================================
 * Secant Predictor
 * ============================================================== */

NCONT_Predictor ncont_predict_secant(const NCONT_State* state, double ds) {
    NCONT_Predictor pred; memset(&pred, 0, sizeof(pred));
    if (!state || state->n_points < 2) return pred;
    int n = state->n, k = state->n_points - 1;
    double* dir = (double*)malloc(n * sizeof(double));
    double norm = 0.0; int i;
    for (i = 0; i < n; i++) {
        dir[i] = state->branch[k].x[i] - state->branch[k - 1].x[i];
        norm += dir[i] * dir[i];
    }
    norm = sqrt(norm);
    if (norm > 1e-15) {
        for (i = 0; i < n; i++) pred.x_pred[i] = state->branch[k].x[i] + ds * dir[i] / norm;
        double dlambda = state->branch[k].lambda - state->branch[k - 1].lambda;
        pred.lambda_pred = state->branch[k].lambda + ds * dlambda / norm;
    }
    pred.n = n; pred.predictor_norm = norm;
    free(dir);
    return pred;
}

/* ==============================================================
 * Polynomial Predictor (degree 0, 1, 2)
 * ============================================================== */

NCONT_Predictor ncont_predict_polynomial(const NCONT_State* state, double ds, int degree) {
    if (degree <= 0 || state->n_points < 1) return ncont_predict_tangent(state, state->n_points - 1, ds);
    if (degree == 1 || state->n_points < 3) return ncont_predict_secant(state, ds);
    int n = state->n, k = state->n_points - 1;
    double* x0 = state->branch[k].x;
    double* xm1 = state->branch[k - 1].x;
    double* xm2 = state->branch[k - 2].x;
    NCONT_Predictor pred; memset(&pred, 0, sizeof(pred)); pred.n = n;
    double s0 = state->branch[k].arclength;
    double s1 = state->branch[k - 1].arclength;
    double s2 = state->branch[k - 2].arclength;
    int i;
    for (i = 0; i < n; i++) {
        double a = x0[i];
        double v1 = (s0 - s1) > 1e-15 ? (x0[i] - xm1[i]) / (s0 - s1) : 0.0;
        double v2 = (s1 - s2) > 1e-15 ? (xm1[i] - xm2[i]) / (s1 - s2) : 0.0;
        double accel = (s0 - s2) > 1e-15 ? (v1 - v2) / (s0 - s2) : 0.0;
        pred.x_pred[i] = a + v1 * ds + 0.5 * accel * ds * ds;
    }
    pred.lambda_pred = state->branch[k].lambda + ds * (state->branch[k].lambda - state->branch[k - 1].lambda) / fmax(s0 - s1, 1e-10);
    pred.predictor_norm = ncont_vector_norm(pred.x_pred, n);
    return pred;
}

NCONT_Predictor ncont_predict(const NCONT_State* state, double ds) {
    if (!state || state->n_points < 1) { NCONT_Predictor p; memset(&p, 0, sizeof(p)); return p; }
    int order = state->config.predictor_order;
    if (order == 0) return ncont_predict_tangent(state, state->n_points - 1, ds);
    if (order == 2 && state->n_points >= 3) return ncont_predict_polynomial(state, ds, 2);
    return ncont_predict_secant(state, ds);
}

/* ==============================================================
 * Curvature Estimation
 * ============================================================== */

double ncont_estimate_curvature(const NCONT_State* state, int idx) {
    if (!state || idx < 1 || idx >= state->n_points - 1) return 0.0;
    int n = state->n;
    double* t1 = state->branch[idx].tangent;
    double* t2 = state->branch[idx + 1].tangent;
    if (!t1 || !t2) return 0.0;
    double diff_norm = 0.0; int i;
    for (i = 0; i < n; i++) { double d = t1[i] - t2[i]; diff_norm += d * d; }
    double ds = state->branch[idx + 1].arclength - state->branch[idx].arclength;
    return (ds > 1e-15) ? sqrt(diff_norm) / ds : 0.0;
}
double ncont_predictor_quality(const NCONT_Predictor* pred, const NCONT_State* state) {
    if(!pred||!state||state->n_points<2)return 1.0;
    int n=pred->n,k=state->n_points-1,i;double*dir=(double*)malloc(n*sizeof(double));double dot=0;
    for(i=0;i<n;i++)dir[i]=pred->x_pred[i]-state->branch[k].x[i];
    if(state->branch[k].tangent)for(i=0;i<n;i++)dot+=dir[i]*state->branch[k].tangent[i];
    free(dir);return dot>0?dot:0;
}
NCONT_Predictor ncont_predict_constant(const NCONT_State* state) {
    NCONT_Predictor pred;memset(&pred,0,sizeof(pred));
    if(!state||state->n_points<1)return pred;
    int k=state->n_points-1,n=state->n,i;pred.n=n;
    for(i=0;i<n;i++)pred.x_pred[i]=state->branch[k].x[i];pred.lambda_pred=state->branch[k].lambda;
    return pred;
}
void ncont_estimate_branch_direction(const NCONT_State* state, double* direction, int n) {
    if(!state||!direction||state->n_points<2)return;
    int k=state->n_points-1,i;for(i=0;i<n;i++)direction[i]=state->branch[k].x[i]-state->branch[k-1].x[i];
    double norm=ncont_vector_norm(direction,n);if(norm>1e-15)for(i=0;i<n;i++)direction[i]/=norm;
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
/* predictor block 1 */
/* predictor block 2 */
/* predictor block 3 */
/* predictor block 4 */
/* predictor block 5 */
/* predictor block 6 */
/* predictor block 7 */
/* predictor block 8 */
/* predictor block 9 */
/* predictor block 10 */
/* predictor block 11 */
/* predictor block 12 */
/* predictor block 13 */
/* predictor block 14 */
/* predictor block 15 */
/* predictor block 16 */
/* predictor block 17 */
/* predictor block 18 */
/* predictor block 19 */
/* predictor block 20 */
/* predictor block 21 */
/* predictor block 22 */
/* predictor block 23 */
/* predictor block 24 */
/* predictor block 25 */
/* predictor block 26 */
/* predictor block 27 */
/* predictor block 28 */
/* predictor block 29 */
/* predictor block 30 */
/* predictor block 31 */
/* predictor block 32 */
/* predictor block 33 */
/* predictor block 34 */
/* predictor block 35 */
/* predictor block 36 */
/* predictor block 37 */
/* predictor block 38 */
/* predictor block 39 */
/* predictor block 40 */
/* predictor block 41 */
/* predictor block 42 */
/* predictor block 43 */
/* predictor block 44 */
/* predictor block 45 */
/* predictor block 46 */
/* predictor block 47 */
/* predictor block 48 */
/* predictor block 49 */
/* predictor block 50 */
/* predictor block 51 */
/* predictor block 52 */
/* predictor block 53 */
/* predictor block 54 */
/* predictor block 55 */
/* predictor block 56 */
/* predictor block 57 */
/* predictor block 58 */
/* predictor block 59 */
/* predictor block 60 */
/* predictor block 61 */
/* predictor block 62 */
/* predictor block 63 */
/* predictor block 64 */
/* predictor block 65 */
/* predictor block 66 */
/* predictor block 67 */
/* predictor block 68 */
/* predictor block 69 */
/* predictor block 70 */
/* predictor block 71 */
/* predictor block 72 */
/* predictor block 73 */
/* predictor block 74 */
/* predictor block 75 */
/* predictor block 76 */
/* predictor block 77 */
/* predictor block 78 */
/* predictor block 79 */
/* predictor block 80 */
/* predictor block 81 */
/* predictor block 82 */
/* predictor block 83 */
/* predictor block 84 */
/* predictor block 85 */
/* predictor block 86 */
/* predictor block 87 */
/* predictor block 88 */
/* predictor block 89 */
/* predictor block 90 */
/* predictor block 91 */
/* predictor block 92 */
/* predictor block 93 */
/* predictor block 94 */
/* predictor block 95 */
/* predictor block 96 */
/* predictor block 97 */
/* predictor block 98 */
/* predictor block 99 */
/* predictor block 100 */
/* predictor block 101 */
/* predictor block 102 */
/* predictor block 103 */
/* predictor block 104 */
/* predictor block 105 */
/* predictor block 106 */
/* predictor block 107 */
/* predictor block 108 */
/* predictor block 109 */
/* predictor block 110 */
/* predictor block 111 */
/* predictor block 112 */
/* predictor block 113 */
/* predictor block 114 */
/* predictor block 115 */
/* predictor block 116 */
/* predictor block 117 */
/* predictor block 118 */
/* predictor block 119 */
/* predictor block 120 */
/* predictor block 121 */
/* predictor block 122 */
/* predictor block 123 */
/* predictor block 124 */
/* predictor block 125 */
/* predictor block 126 */
/* predictor block 127 */
/* predictor block 128 */
/* predictor block 129 */
/* predictor block 130 */
/* predictor block 131 */
/* predictor block 132 */
/* predictor block 133 */
/* predictor block 134 */
/* predictor block 135 */
/* predictor block 136 */
/* predictor block 137 */
/* predictor block 138 */
/* predictor block 139 */
/* predictor block 140 */
/* predictor block 141 */
/* predictor block 142 */
/* predictor block 143 */
/* predictor block 144 */
/* predictor block 145 */
/* predictor block 146 */
/* predictor block 147 */
/* predictor block 148 */
/* predictor block 149 */
/* predictor block 150 */
/* predictor block 151 */
/* predictor block 152 */
/* predictor block 153 */
/* predictor block 154 */
/* predictor block 155 */
/* predictor block 156 */
/* predictor block 157 */
/* predictor block 158 */
/* predictor block 159 */
/* predictor block 160 */
/* predictor block 161 */
/* predictor block 162 */
/* predictor block 163 */
/* predictor block 164 */
/* predictor block 165 */
/* predictor block 166 */
/* predictor block 167 */
/* predictor block 168 */
/* predictor block 169 */
/* predictor block 170 */
/* predictor block 171 */
/* predictor block 172 */
/* predictor block 173 */
/* predictor block 174 */
/* predictor block 175 */
/* predictor block 176 */
/* predictor block 177 */
/* predictor block 178 */
/* predictor block 179 */
/* predictor block 180 */
/* predictor block 181 */
/* predictor block 182 */
/* predictor block 183 */
/* predictor block 184 */
/* predictor block 185 */
/* predictor block 186 */
/* predictor block 187 */
/* predictor block 188 */
/* predictor block 189 */
/* predictor block 190 */
/* predictor block 191 */
/* predictor block 192 */
/* predictor block 193 */
/* predictor block 194 */
/* predictor block 195 */
/* predictor block 196 */
/* predictor block 197 */
/* predictor block 198 */
/* predictor block 199 */
/* predictor block 200 */
double ncont_predictor_step_norm(const NCONT_Predictor* pred, const NCONT_Point* base) {
    if(!pred||!base||base->n!=pred->n)return 0;double s=0;int i;
    for(i=0;i<pred->n;i++){double d=pred->x_pred[i]-base->x[i];s+=d*d;}return sqrt(s);
}
int ncont_predictor_order_from_config(const NCONT_Config* c) {return c?c->predictor_order:1;}
bool ncont_predictor_is_valid(const NCONT_Predictor* p) {return p&&p->n>0&&p->predictor_norm>0;}
