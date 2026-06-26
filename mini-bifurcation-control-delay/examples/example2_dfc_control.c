/* Delayed feedback control: compute optimal Pyragas gain for Rossler UPO */
#include "../include/bcd_control.h"
#include "../include/bcd_core.h"
#include <stdio.h>
#include <math.h>

static void rossler_rhs(int n, const double* x, const double* xd, double t, double* dx, void* p) {
    (void)xd; (void)t; double* pars = (double*)p;
    double a = pars[0], b = pars[1], c = pars[2];
    dx[0] = -x[1] - x[2];
    dx[1] = x[0] + a * x[1];
    dx[2] = b + x[2] * (x[0] - c);
}

int main(void) {
    printf("=== Delayed Feedback Control: Optimal Gain Design ===\n\n");

    double pars[] = {0.2, 0.2, 5.7};
    double tau_upo = 5.88;  /* known UPO period */
    double x0[] = {1.0, 0.5, 2.0};

    /* Setup DDE system */
    BCDDDE* dde = bcd_dde_create(3, tau_upo, 600);
    bcd_dde_set_history(dde, x0, 3);

    /* Compute optimal Pyragas gain */
    double K_opt = bcd_pyragas_optimal_gain(dde, rossler_rhs, pars, 6.283185307 / tau_upo);
    printf("Optimal Pyragas gain K = %.6f (omega=%.4f)\n", K_opt, 6.283185307 / tau_upo);

    /* Design full control */
    BCDControlConfig cfg = bcd_control_config_default();
    cfg.type = BCD_CTRL_PYRAGAS;
    cfg.gain = K_opt;
    cfg.delay = tau_upo;
    cfg.dim = 3;
    BCDControlResult* cr = bcd_control_design(dde, rossler_rhs, pars, &cfg);

    if (cr) {
        printf("Control energy: %.6f, settling time: %.4f\n",
               cr->energy, cr->settling_time);
        printf("Steady-state error: %.6f, stabilized: %s\n",
               cr->steady_error, bcd_control_is_stabilized(cr, 0.01) ? "YES" : "NO");
        double eff = bcd_control_energy_efficiency(cr);
        printf("Energy efficiency: %.4f\n", eff);

        /* Compare with extended DFC */
        double K_ext = bcd_extended_dfc_gain(dde, rossler_rhs, pars, 3);
        printf("Extended DFC (3 terms) gain: %.6f\n", K_ext);
        bcd_control_result_free(cr);
    }

    /* Nyquist analysis */
    double A_test[] = {0.0, 0.0, 0.0, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0};
    double B_test[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    double nyq = bcd_nyquist_criterion_delay(A_test, B_test, 3, tau_upo);
    printf("Nyquist criterion (delay): %.4f\n", nyq);

    bcd_dde_free(dde);
    printf("\nExample 2 PASSED\n");
    return 0;
}
