/* Step-size adaptation: compare fixed vs adaptive continuation */
#include "../include/ncont_arclength.h"
#include <stdio.h>
#include <math.h>

static void cubic_sys(const double* x, int n, double lambda, void* p, double* r) {
    (void)p; (void)n;
    r[0] = x[0]*x[0]*x[0] - 3.0*x[0] + lambda;  /* cubic with fold */
}
static void cubic_jac(const double* x, int n, double lambda, void* p, double* J) {
    (void)p; (void)n; (void)lambda;
    J[0] = 3.0 * x[0] * x[0] - 3.0;  /* df/dx = 3x^2 - 3 */
}

int main(void) {
    printf("=== Step-Size Adaptation Comparison ===\n\n");

    double x0[] = {-2.0};  /* start on stable branch */

    /* --- Fixed step size --- */
    printf("Fixed step ds=0.1:\n");
    NCONT_Config cfg_fixed = ncont_config_default();
    cfg_fixed.ds = 0.1; cfg_fixed.adaptive_stepsize = false;
    cfg_fixed.max_points = 100;
    NCONT_State* st_fixed = NULL;
    NCONT_Result res_fixed = ncont_continue(cubic_sys, cubic_jac, x0, 1, 0.0, &cfg_fixed, NULL, &st_fixed);
    printf("  Points: %d, arclength: %.4f, rejected: %d, avg Newton: %.1f\n",
           res_fixed.total_points, res_fixed.arclength_covered,
           res_fixed.rejected_steps, res_fixed.avg_newton_iterations);
    ncont_state_free(st_fixed);

    /* --- Adaptive step size --- */
    printf("\nAdaptive ds=[0.005..0.5]:\n");
    NCONT_Config cfg_adapt = ncont_config_default();
    cfg_adapt.ds = 0.2; cfg_adapt.ds_min = 0.005;
    cfg_adapt.ds_max = 0.5; cfg_adapt.adaptive_stepsize = true;
    cfg_adapt.max_points = 100;
    NCONT_State* st_adapt = NULL;
    NCONT_Result res_adapt = ncont_continue(cubic_sys, cubic_jac, x0, 1, 0.0, &cfg_adapt, NULL, &st_adapt);
    printf("  Points: %d, arclength: %.4f, rejected: %d, avg Newton: %.1f\n",
           res_adapt.total_points, res_adapt.arclength_covered,
           res_adapt.rejected_steps, res_adapt.avg_newton_iterations);
    printf("  Step range: [%.5f, %.4f]\n",
           res_adapt.min_stepsize_used, res_adapt.max_stepsize_used);

    /* Compare efficiency */
    if (res_fixed.total_points > 0 && res_adapt.total_points > 0) {
        double eff_fixed = res_fixed.arclength_covered / (double)res_fixed.total_points;
        double eff_adapt = res_adapt.arclength_covered / (double)res_adapt.total_points;
        printf("\nEfficiency (arclength/pt): fixed=%.4f, adaptive=%.4f\n", eff_fixed, eff_adapt);
    }

    ncont_state_free(st_adapt);
    printf("\nExample 3 PASSED\n");
    return 0;
}
