/* bcd_app2.c -- L7-2: Laser stabilization via Pyragas DFC (Bielawski 1994).
 * A CO2 laser with feedback exhibits chaotic intensity fluctuations.
 * Pyragas delayed feedback u(t) = K*(I(t-tau) - I(t)), with tau matched
 * to the unstable periodic orbit period, stabilizes the laser output.
 * This application computes the optimal control parameters. */
#include "../include/bcd_control.h"
#include "../include/bcd_core.h"
#include <math.h>
#include <stdio.h>

static void laser_dde_rhs(int n, const double* x, const double* xd, double t, double* dx, void* p) {
    (void)t; (void)xd;
    double* pars = (double*)p;
    double I = x[0], N = x[1];
    double kappa = pars[0], gamma = pars[1], A = pars[2];
    dx[0] = -kappa * I + A * I * N;  /* intensity */
    dx[1] = -gamma * N + 1.0 - I * N;  /* population inversion */
}

int bcd_laser_dfc(double kappa, double gamma, double A, double* opt_gain,
                  double* opt_delay, double* stabilization_quality) {
    double tau_guess = 2.0 * M_PI / sqrt(kappa * gamma * A);  /* relaxation osc period */
    double pars[] = {kappa, gamma, A};

    BCDDDE* dde = bcd_dde_create(2, tau_guess, 500);
    double x0[] = {0.5, 0.5};
    bcd_dde_set_history(dde, x0, 2);

    /* Compute optimal Pyragas gain */
    *opt_gain = bcd_pyragas_optimal_gain(dde, laser_dde_rhs, pars, 2.0*M_PI/tau_guess);
    *opt_delay = tau_guess;

    /* Design full DFC control */
    BCDControlConfig cfg = bcd_control_config_default();
    cfg.type = BCD_CTRL_PYRAGAS;
    cfg.gain = *opt_gain;
    cfg.delay = tau_guess;
    cfg.dim = 2;
    BCDControlResult* cr = bcd_control_design(dde, laser_dde_rhs, pars, &cfg);

    if (cr) {
        *stabilization_quality = cr->success ? (1.0 / (1.0 + cr->steady_error)) : 0.0;
        printf("[bcd_app2] Laser DFC stabilization:\n");
        printf("  kappa=%.2f gamma=%.2f A=%.2f\n", kappa, gamma, A);
        printf("  Optimal: K=%.6f, tau=%.6f s\n", *opt_gain, *opt_delay);
        printf("  Settling time: %.4f, steady error: %.6f\n", cr->settling_time, cr->steady_error);
        printf("  Stabilized: %s (quality=%.4f)\n",
               bcd_control_is_stabilized(cr, 0.01) ? "YES" : "NO", *stabilization_quality);

        /* Compare with extended DFC */
        double K_ext = bcd_extended_dfc_gain(dde, laser_dde_rhs, pars, 5);
        double eff = bcd_control_energy_efficiency(cr);
        printf("  Extended DFC (5 terms): K=%.4f, efficiency=%.4f\n", K_ext, eff);
        bcd_control_result_free(cr);
    }

    bcd_dde_free(dde);
    return 0;
}

int main(void) {
    double K, tau, quality;
    bcd_laser_dfc(10.0, 1.0, 0.8, &K, &tau, &quality);
    printf("=== bcd_app2: Laser DFC Stabilization ===\n");
    printf("K_opt=%.4f tau=%.4f quality=%.4f\n", K, tau, quality);
    return 0;
}
