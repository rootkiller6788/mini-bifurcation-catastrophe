/* bcd_app1.c -- L7: Machining chatter suppression via DFC (Tobias 1965).
 * Machine tool regenerative chatter is a DDE: mxddot + c*xdot + k*x = -Kc*(x(t)-x(t-tau))
 * where tau = 2*pi/Omega (spindle period). By applying delayed feedback
 * control with appropriate gain and delay matching, chatter can be
 * suppressed — a practical application of bifurcation control. */
#include "../include/bcd_control.h"
#include "../include/bcd_core.h"
#include <math.h>
#include <stdio.h>

static void chatter_dde(int n, const double* x, const double* xd, double t, double* dx, void* p) {
    (void)t;
    double* pars = (double*)p;
    double m = pars[0], c = pars[1], k = pars[2], Kc = pars[3];
    double x_undelayed = x[0];
    double x_delayed = xd[0];
    dx[0] = x[1];  /* dx/dt = v */
    dx[1] = (-c * x[1] - k * x_undelayed - Kc * (x_undelayed - x_delayed)) / m;
}

int bcd_chatter_control(double spindle_speed_rpm, double chatter_freq,
                        double m, double c, double k, double* optimal_gain,
                        double* stability_lobe) {
    double omega = spindle_speed_rpm * 2.0 * M_PI / 60.0;
    double tau = 2.0 * M_PI / omega;  /* one revolution delay */
    double pars[] = {m, c, k, 0.0};   /* Kc will be determined */

    BCDDDE* dde = bcd_dde_create(2, tau, 200);
    double x0[] = {0.01, 0.0};
    bcd_dde_set_history(dde, x0, 2);

    /* Compute optimal Pyragas gain for chatter frequency */
    double omega_c = chatter_freq * 2.0 * M_PI;
    *optimal_gain = bcd_pyragas_optimal_gain(dde, chatter_dde, pars, omega_c);

    /* Stability lobe: minimum chip width before chatter */
    double k_crit = c / (2.0 * m * omega_c);
    *stability_lobe = k_crit / fabs(*optimal_gain);

    printf("[bcd_app1] Machining chatter control:\n");
    printf("  Spindle: %.0f RPM (tau=%.4f s), chatter freq: %.1f Hz\n",
           spindle_speed_rpm, tau, chatter_freq);
    printf("  Optimal DFC gain K=%.4f (negative => feedback stabilization)\n", *optimal_gain);
    printf("  Stability lobe: min depth=%.4f mm\n", *stability_lobe);

    /* Simulate with and without control */
    for (int ctrl_on = 0; ctrl_on <= 1; ctrl_on++) {
        BCDDDE* sim = bcd_dde_create(2, tau, 500);
        bcd_dde_set_history(sim, x0, 2);
        double Kc = ctrl_on ? *optimal_gain : 0.0;
        double sim_pars[] = {m, c, k, Kc};
        for (int s = 0; s < 400; s++)
            bcd_dde_step(sim, chatter_dde, sim_pars, 0.01);
        double xf = sim->buffer[sim->buf_pos];
        printf("  Control %s: final |x|=%.6f %s\n",
               ctrl_on ? "ON " : "OFF", fabs(xf),
               fabs(xf) < 0.001 ? "(stable)" : "(chatter)");
        bcd_dde_free(sim);
    }

    bcd_dde_free(dde);
    return 0;
}

int main(void) {
    double K_opt, lobe;
    bcd_chatter_control(3000.0, 120.0, 0.5, 2.0, 5000.0, &K_opt, &lobe);
    printf("=== bcd_app1: Chatter Suppression ===\n");
    printf("K_opt=%.4f lobe=%.4f\n", K_opt, lobe);
    return 0;
}
