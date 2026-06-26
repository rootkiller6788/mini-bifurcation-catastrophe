/* Numerical continuation: trace fold bifurcation via arc-length method */
#include "../include/ncont_arclength.h"
#include "../include/ncont_bifurcation.h"
#include <stdio.h>
#include <math.h>

static void fold_sys(const double* x, int n, double lambda, void* p, double* r) {
    (void)p; (void)n;
    r[0] = x[0] * x[0] + lambda;  /* f(x,lambda) = x^2 + lambda = 0 */
}
static void fold_jac(const double* x, int n, double lambda, void* p, double* J) {
    (void)p; (void)n; (void)lambda;
    J[0] = 2.0 * x[0];  /* df/dx = 2x */
}

int main(void) {
    printf("=== Numerical Continuation: Fold Bifurcation ===\n\n");

    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.1;
    cfg.ds_min = 0.001;
    cfg.ds_max = 0.5;
    cfg.newton_tol = 1e-8;
    cfg.newton_max_iter = 30;
    cfg.max_points = 80;
    cfg.detect_bifurcations = true;
    cfg.adaptive_stepsize = true;

    double x0[] = {1.0};
    NCONT_State* state = NULL;

    printf("Starting at x=%.2f, lambda=%.2f\n", x0[0], -1.0);
    NCONT_Result res = ncont_continue(fold_sys, fold_jac, x0, 1, -1.0, &cfg, NULL, &state);

    printf("\nResult: %d points, arclength=%.4f, exit=%d\n",
           res.total_points, res.arclength_covered, res.exit_flag);
    printf("Accepted: %d, rejected: %d, bif detected: %d\n",
           res.accepted_points, res.rejected_steps, res.n_bifurcations_detected);
    printf("Step size: min=%.6f, max=%.4f, avg Newton iters=%.2f\n",
           res.min_stepsize_used, res.max_stepsize_used, res.avg_newton_iterations);

    if (state && state->n_points > 0) {
        printf("\nBranch trace (first/last/mid):\n");
        for (int i = 0; i < state->n_points; i += (state->n_points > 10 ? state->n_points/5 : 1)) {
            printf("  pt[%d]: x=%.4f lambda=%.4f\n", i,
                   state->branch[i].x[0], state->branch[i].lambda);
        }
    }

    /* Bifurcation scan */
    if (state && res.n_bifurcations_detected > 0) {
        NCONT_Bifurcation bif_list[10];
        int nb = ncont_scan_bifurcations(state, fold_jac, 1, NULL, bif_list, 10);
        printf("\nBifurcations found: %d\n", nb);
        for (int i = 0; i < nb; i++) {
            printf("  %s at lambda=%.4f, test=%.6f\n",
                   ncont_bifurcation_type_name(bif_list[i].type),
                   bif_list[i].lambda, bif_list[i].test_value);
        }
    }

    ncont_state_free(state);
    printf("\nExample 1 PASSED\n");
    return 0;
}
