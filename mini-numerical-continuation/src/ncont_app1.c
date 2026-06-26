/* ncont_app1.c -- L7: Aircraft wing load continuation (Crisfield 1981, arc-length).
 * Wing structures exhibit snap-through buckling under aerodynamic loading.
 * The load-displacement curve has a limit point (fold) that cannot be traced
 * by simple load incrementation. Arc-length continuation passes through the
 * limit point to reveal the unstable branch and post-buckling behavior. */
#include "../include/ncont_arclength.h"
#include <math.h>
#include <stdio.h>

/* Wing panel: cubic-softening spring model F(u) = k1*u - k3*u^3 = lambda */
static void wing_sys(const double* x, int n, double lambda, void* p, double* r) {
    double* pars = (double*)p;  /* k1, k3 */
    r[0] = pars[0] * x[0] - pars[1] * x[0] * x[0] * x[0] - lambda;
}
static void wing_jac(const double* x, int n, double lambda, void* p, double* J) {
    (void)lambda; double* pars = (double*)p;
    J[0] = pars[0] - 3.0 * pars[1] * x[0] * x[0];
}

int ncont_wing_limit(double k1, double k3, double* P_crit, double* disp_curve, int n_pts) {
    double pars[] = {k1, k3};
    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.1; cfg.ds_max = 0.3; cfg.ds_min = 0.002;
    cfg.max_points = n_pts; cfg.detect_bifurcations = true;
    cfg.adaptive_stepsize = true;

    double x0[] = {0.5};  /* start on stable branch */
    NCONT_State* state = NULL;

    NCONT_Result res = ncont_continue(wing_sys, wing_jac, x0, 1, 0.0, &cfg, pars, &state);

    printf("[ncont_app1] Wing load continuation: k1=%.2f k3=%.2f\n", k1, k3);
    printf("  Points: %d, arclength: %.4f, bif detected: %d\n",
           res.total_points, res.arclength_covered, res.n_bifurcations_detected);

    /* Find the limit point (max lambda) */
    double lambda_max = -1e10;
    int max_idx = -1;
    for (int i = 0; i < state->n_points; i++) {
        if (state->branch[i].lambda > lambda_max) {
            lambda_max = state->branch[i].lambda;
            max_idx = i;
        }
        if (i < n_pts) disp_curve[i] = state->branch[i].x[0];
    }
    *P_crit = lambda_max;

    printf("  Limit point: lambda_crit=%.6f at u=%.6f (idx %d)\n",
           lambda_max, state->branch[max_idx].x[0], max_idx);

    /* Print load-displacement curve summary */
    printf("  Load-displacement (first/last/mid):\n");
    for (int i = 0; i < state->n_points; i += state->n_points/5) {
        printf("    u=%.4f lambda=%.4f\n", state->branch[i].x[0], state->branch[i].lambda);
    }

    ncont_state_free(state);
    return res.exit_flag == 0 ? 0 : 1;
}

int main(void) {
    double P_crit, disp[100];
    ncont_wing_limit(1.0, 0.25, &P_crit, disp, 80);
    printf("=== ncont_app1: Wing Load Limit Point ===\n");
    printf("P_crit=%.4f\n", P_crit);
    return 0;
}
