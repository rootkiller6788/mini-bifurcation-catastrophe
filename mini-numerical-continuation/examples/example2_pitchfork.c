/* Numerical continuation: pitchfork branch switching */
#include "../include/ncont_branch.h"
#include "../include/ncont_bifurcation.h"
#include "../include/ncont_predictor.h"
#include <stdio.h>
#include <math.h>

static void pitch_sys(const double* x, int n, double lambda, void* p, double* r) {
    (void)p; (void)n;
    r[0] = lambda * x[0] - x[0] * x[0] * x[0];  /* f = lambda*x - x^3 */
}
static void pitch_jac(const double* x, int n, double lambda, void* p, double* J) {
    (void)p; (void)n;
    J[0] = lambda - 3.0 * x[0] * x[0];  /* df/dx = lambda - 3x^2 */
}

int main(void) {
    printf("=== Continuation: Pitchfork Branch Switching ===\n\n");

    /* First: trace trivial branch x=0 */
    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.1; cfg.max_points = 60;
    cfg.detect_bifurcations = true;
    cfg.adaptive_stepsize = true;
    cfg.newton_tol = 1e-8;

    double x0[] = {0.5};
    NCONT_State* state = NULL;
    printf("Tracing trivial branch from x=0.5, lambda=-2.0...\n");
    NCONT_Result res = ncont_continue(pitch_sys, pitch_jac, x0, 1, -2.0, &cfg, NULL, &state);
    printf("Branch: %d points, exit=%d\n", res.total_points, res.exit_flag);

    /* Detect bifurcation */
    NCONT_Bifurcation bif_list[5];
    int nb = ncont_scan_bifurcations(state, pitch_jac, 1, NULL, bif_list, 5);
    printf("Bifurcations found: %d\n", nb);
    for (int i = 0; i < nb; i++) {
        printf("  %s at lambda=%.4f (point idx %d)\n",
               ncont_bifurcation_type_name(bif_list[i].type),
               bif_list[i].lambda, bif_list[i].point_index);
    }

    /* Branch switching at detected bifurcation */
    if (nb > 0) {
        NCONT_Point* bif_pt = &state->branch[bif_list[0].point_index];
        printf("\nSwitching branch at lambda=%.4f...\n", bif_pt->lambda);
        NCONT_BranchSwitch sw = ncont_switch_branch(pitch_sys, pitch_jac,
            &bif_list[0], bif_pt, 1, NULL, 0.05);
        printf("Switch: success=%s, dir=%d, lambda_new=%.4f\n",
               sw.success ? "YES" : "NO", sw.direction, sw.lambda_new);

        /* Trace new branch */
        if (sw.success) {
            NCONT_State* state2 = NULL;
            NCONT_Config cfg2 = ncont_config_default();
            cfg2.ds = 0.05; cfg2.max_points = 40;
            cfg2.arclength_max = 3.0;
            NCONT_Result res2 = ncont_continue(pitch_sys, pitch_jac, sw.x_new, 1,
                sw.lambda_new, &cfg2, NULL, &state2);
            printf("New branch: %d points, arclength=%.4f\n",
                   res2.total_points, res2.arclength_covered);
            ncont_state_free(state2);
        }
    }

    ncont_state_free(state);
    printf("\nExample 2 PASSED\n");
    return 0;
}
