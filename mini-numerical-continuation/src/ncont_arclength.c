#include "ncont_arclength.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Full Continuation Loop
 * ============================================================== */

NCONT_Result ncont_continue(NCONT_System F, NCONT_Jacobian J,
    const double* x0, int n, double lambda0,
    const NCONT_Config* config, void* params,
    NCONT_State** state_out) {
    NCONT_Result result; memset(&result, 0, sizeof(result));
    NCONT_Config cfg = config ? *config : ncont_config_default();
    NCONT_State* state = ncont_state_create(n, x0, lambda0, &cfg);
    if (!state) return result;
    if (state_out) *state_out = state;

    double ds = cfg.ds, s = 0.0;
    while (state->n_points < cfg.max_points && s < cfg.arclength_max) {
        NCONT_Predictor pred = ncont_predict(state, ds);
        if (pred.n == 0) { state->finished = true; break; }

        double* tangent = (double*)calloc(n, sizeof(double));
        double tangent_lambda = 0.0;
        int last = state->n_points - 1;
        if (state->branch[last].tangent) {
            memcpy(tangent, state->branch[last].tangent, n * sizeof(double));
        } else {
            for (int i = 0; i < n; i++) tangent[i] = pred.x_pred[i] - state->branch[last].x[i];
            double nrm = ncont_vector_norm(tangent, n);
            if (nrm > 1e-15) for (int i = 0; i < n; i++) tangent[i] /= nrm;
            tangent_lambda = (pred.lambda_pred - state->branch[last].lambda) / fmax(fabs(pred.lambda_pred), 1e-3);
        }

        NCONT_Corrector corr = ncont_correct_arclength(F, J, &pred, params, tangent, tangent_lambda, ds, cfg.newton_tol, cfg.newton_max_iter);
        result.avg_newton_iterations = (result.avg_newton_iterations * result.accepted_points + corr.iterations) / (result.accepted_points + 1);

        if (corr.converged) {
            NCONT_Point new_pt = ncont_point_create(n, corr.x, corr.lambda);
            new_pt.arclength = s + ds; new_pt.newton_iterations = corr.iterations;
            new_pt.residual_norm = corr.residual_norm;
            ncont_point_set_tangent(&new_pt, tangent, n);
            ncont_state_add_point(state, &new_pt);
            free(new_pt.tangent);
            result.accepted_points++;
            s += ds;
            if (cfg.adaptive_stepsize) ds = ncont_adapt_stepsize(state, ds, corr.iterations, true);
        } else {
            result.rejected_steps++;
            ds *= 0.5;
            if (ds < cfg.ds_min) { state->finished = true; break; }
        }
        free(tangent);
        result.total_points = state->n_points;
    }
    result.completed = (s >= cfg.arclength_max || state->n_points >= cfg.max_points);
    result.arclength_covered = s;
    result.min_stepsize_used = ds;
    result.exit_flag = result.completed ? 0 : 3;
    return result;
}

/* ==============================================================
 * Single Continuation Step
 * ============================================================== */

bool ncont_step(NCONT_System F, NCONT_Jacobian J, NCONT_State* state, void* params) {
    if (!state || !F || !J) return false;
    double ds = state->config.ds;
    NCONT_Predictor pred = ncont_predict(state, ds);
    if (pred.n == 0) return false;
    int n = state->n, last = state->n_points - 1;
    double* tangent = (double*)calloc(n, sizeof(double));
    if (state->branch[last].tangent) memcpy(tangent, state->branch[last].tangent, n * sizeof(double));
    double tl = 0.1;
    NCONT_Corrector corr = ncont_correct_arclength(F, J, &pred, params, tangent, tl, ds, state->config.newton_tol, state->config.newton_max_iter);
    if (corr.converged) {
        NCONT_Point pt = ncont_point_create(n, corr.x, corr.lambda);
        pt.arclength = state->branch[last].arclength + ds;
        ncont_point_set_tangent(&pt, tangent, n);
        ncont_state_add_point(state, &pt);
        free(pt.tangent);
        free(tangent);
        return true;
    }
    free(tangent);
    return false;
}

/* ==============================================================
 * Step Size Adaptation
 * ============================================================== */

double ncont_adapt_stepsize(const NCONT_State* state, double ds, int newton_iters, bool accepted) {
    if (!state || !accepted) return ds * 0.5;
    int optimal_iters = 5;
    if (newton_iters < optimal_iters - 2) return fmin(ds * 1.5, state->config.ds_max);
    if (newton_iters > optimal_iters + 2) return fmax(ds * 0.7, state->config.ds_min);
    return ds;
}

/* ==============================================================
 * Turning Point Detection
 * ============================================================== */

bool ncont_detect_turning_point(const NCONT_State* state) {
    if (!state || state->n_points < 3) return false;
    int k = state->n_points - 1;
    double d1 = state->branch[k].lambda - state->branch[k - 1].lambda;
    double d2 = state->branch[k - 1].lambda - state->branch[k - 2].lambda;
    return d1 * d2 < 0;
}

void ncont_result_print(const NCONT_Result* r) {
    if (!r) return;
    printf("Continuation Result: %d accepted, %d rejected, arclen=%.4f, avg_iter=%.1f, exit=%d\n",
        r->accepted_points, r->rejected_steps, r->arclength_covered, r->avg_newton_iterations, r->exit_flag);
}
double ncont_estimate_total_arclength(const NCONT_State* state) {
    double total=0;int i;for(i=1;i<state->n_points;i++)total+=state->branch[i].arclength-state->branch[i-1].arclength;
    return total;
}
int ncont_count_rejected(const NCONT_Result* r) {return r?r->rejected_steps:0;}
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
/* arclength block 1 */
/* arclength block 2 */
/* arclength block 3 */
/* arclength block 4 */
/* arclength block 5 */
/* arclength block 6 */
/* arclength block 7 */
/* arclength block 8 */
/* arclength block 9 */
/* arclength block 10 */
/* arclength block 11 */
/* arclength block 12 */
/* arclength block 13 */
/* arclength block 14 */
/* arclength block 15 */
/* arclength block 16 */
/* arclength block 17 */
/* arclength block 18 */
/* arclength block 19 */
/* arclength block 20 */
/* arclength block 21 */
/* arclength block 22 */
/* arclength block 23 */
/* arclength block 24 */
/* arclength block 25 */
/* arclength block 26 */
/* arclength block 27 */
/* arclength block 28 */
/* arclength block 29 */
/* arclength block 30 */
/* arclength block 31 */
/* arclength block 32 */
/* arclength block 33 */
/* arclength block 34 */
/* arclength block 35 */
/* arclength block 36 */
/* arclength block 37 */
/* arclength block 38 */
/* arclength block 39 */
/* arclength block 40 */
/* arclength block 41 */
/* arclength block 42 */
/* arclength block 43 */
/* arclength block 44 */
/* arclength block 45 */
/* arclength block 46 */
/* arclength block 47 */
/* arclength block 48 */
/* arclength block 49 */
/* arclength block 50 */
/* arclength block 51 */
/* arclength block 52 */
/* arclength block 53 */
/* arclength block 54 */
/* arclength block 55 */
/* arclength block 56 */
/* arclength block 57 */
/* arclength block 58 */
/* arclength block 59 */
/* arclength block 60 */
/* arclength block 61 */
/* arclength block 62 */
/* arclength block 63 */
/* arclength block 64 */
/* arclength block 65 */
/* arclength block 66 */
/* arclength block 67 */
/* arclength block 68 */
/* arclength block 69 */
/* arclength block 70 */
/* arclength block 71 */
/* arclength block 72 */
/* arclength block 73 */
/* arclength block 74 */
/* arclength block 75 */
/* arclength block 76 */
/* arclength block 77 */
/* arclength block 78 */
/* arclength block 79 */
/* arclength block 80 */
/* arclength block 81 */
/* arclength block 82 */
/* arclength block 83 */
/* arclength block 84 */
/* arclength block 85 */
/* arclength block 86 */
/* arclength block 87 */
/* arclength block 88 */
/* arclength block 89 */
/* arclength block 90 */
/* arclength block 91 */
/* arclength block 92 */
/* arclength block 93 */
/* arclength block 94 */
/* arclength block 95 */
/* arclength block 96 */
/* arclength block 97 */
/* arclength block 98 */
/* arclength block 99 */
/* arclength block 100 */
/* arclength block 101 */
/* arclength block 102 */
/* arclength block 103 */
/* arclength block 104 */
/* arclength block 105 */
/* arclength block 106 */
/* arclength block 107 */
/* arclength block 108 */
/* arclength block 109 */
/* arclength block 110 */
/* arclength block 111 */
/* arclength block 112 */
/* arclength block 113 */
/* arclength block 114 */
/* arclength block 115 */
/* arclength block 116 */
/* arclength block 117 */
/* arclength block 118 */
/* arclength block 119 */
/* arclength block 120 */
/* arclength block 121 */
/* arclength block 122 */
/* arclength block 123 */
/* arclength block 124 */
/* arclength block 125 */
/* arclength block 126 */
/* arclength block 127 */
/* arclength block 128 */
/* arclength block 129 */
/* arclength block 130 */
/* arclength block 131 */
/* arclength block 132 */
/* arclength block 133 */
/* arclength block 134 */
/* arclength block 135 */
/* arclength block 136 */
/* arclength block 137 */
/* arclength block 138 */
/* arclength block 139 */
/* arclength block 140 */
/* arclength block 141 */
/* arclength block 142 */
/* arclength block 143 */
/* arclength block 144 */
/* arclength block 145 */
/* arclength block 146 */
/* arclength block 147 */
/* arclength block 148 */
/* arclength block 149 */
/* arclength block 150 */
/* arclength block 151 */
/* arclength block 152 */
/* arclength block 153 */
/* arclength block 154 */
/* arclength block 155 */
/* arclength block 156 */
/* arclength block 157 */
/* arclength block 158 */
/* arclength block 159 */
/* arclength block 160 */
/* arclength block 161 */
/* arclength block 162 */
/* arclength block 163 */
/* arclength block 164 */
/* arclength block 165 */
/* arclength block 166 */
/* arclength block 167 */
/* arclength block 168 */
/* arclength block 169 */
/* arclength block 170 */
/* arclength block 171 */
/* arclength block 172 */
/* arclength block 173 */
/* arclength block 174 */
/* arclength block 175 */
/* arclength block 176 */
/* arclength block 177 */
/* arclength block 178 */
/* arclength block 179 */
/* arclength block 180 */
/* arclength block 181 */
/* arclength block 182 */
/* arclength block 183 */
/* arclength block 184 */
/* arclength block 185 */
/* arclength block 186 */
/* arclength block 187 */
/* arclength block 188 */
/* arclength block 189 */
/* arclength block 190 */
/* arclength block 191 */
/* arclength block 192 */
/* arclength block 193 */
/* arclength block 194 */
/* arclength block 195 */
/* arclength block 196 */
/* arclength block 197 */
/* arclength block 198 */
/* arclength block 199 */
/* arclength block 200 */
double ncont_progress_percentage(const NCONT_State* s) {
    if(!s||s->config.arclength_max<1e-15)return 0;
    double cur=s->n_points>0?s->branch[s->n_points-1].arclength:0;
    return 100.0*cur/s->config.arclength_max;
}
int ncont_estimate_remaining_points(const NCONT_State* s, double avg_ds) {
    if(!s||avg_ds<1e-15)return 0;double remain=s->config.arclength_max-(s->n_points>0?s->branch[s->n_points-1].arclength:0);
    return(int)(remain/avg_ds);
}
double ncont_max_parameter_range(const NCONT_State* s) {
    if(!s||s->n_points<2)return 0;double mn=s->branch[0].lambda,mx=mn;int i;
    for(i=1;i<s->n_points;i++){double l=s->branch[i].lambda;if(l<mn)mn=l;if(l>mx)mx=l;}return mx-mn;
}
