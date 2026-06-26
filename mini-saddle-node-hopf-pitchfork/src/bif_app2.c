/* bif_app2.c -- L7-2: Laser threshold Hopf (Haken 1975, single-mode laser).
 * Rate equations: dE/dt = -kappa*E + g*E*N, dN/dt = -gamma*N + P - g*E^2*N
 * At threshold P_th = kappa*gamma/g, the off-state (E=0,N=P/gamma) loses
 * stability via a Hopf bifurcation when relaxation oscillation frequency
 * crosses zero damping boundary. */
#include "../include/hopf_bifurcation.h"
#include <math.h>
#include <stdio.h>

static void laser_ode(const double* x, double* params, double* dx, int n) {
    (void)n;
    double kappa = params[0], gamma = params[1], g = params[2], P = params[3];
    dx[0] = -kappa * x[0] + g * x[0] * x[1];      /* field amplitude */
    dx[1] = -gamma * x[1] + P - g * x[0] * x[0] * x[1];  /* population inversion */
}

int bif_laser_threshold(double kappa, double gamma, double g, double P_range[2],
                        double* P_th, double* l1_coeff) {
    double x0[] = {0.05, 0.5};
    double p[] = {kappa, gamma, g, P_range[0]};

    HopfBifurcation* hb = hopf_create();
    int rc = hopf_detect(laser_ode, x0, p, 2, 3, P_range[0], P_range[1],
                          0.02 * (P_range[1] - P_range[0]), 80, 1e-8, hb);

    if (rc >= 0) {
        *P_th = hb->r_critical;
        *l1_coeff = hb->first_lyapunov;

        printf("[bif_app2] Laser threshold: P_th=%.4f (expected %.4f)\n",
               hb->r_critical, kappa * gamma / g);
        printf("  Critical freq omega=%.4f, l1=%.6f, type=%s\n",
               hb->omega_critical, hb->first_lyapunov,
               hb->is_supercritical ? "SUPERCRITICAL" : "SUBCRITICAL");

        /* Compute relaxation oscillation frequency */
        double omega_RO = sqrt(2.0 * kappa * g * hb->r_critical / gamma);
        printf("  Relaxation osc freq: %.4f rad/s\n", omega_RO);

        /* Detect limit cycle above threshold */
        double amp, period;
        double x_lc[] = {0.2, 0.3};
        double p_lc[] = {kappa, gamma, g, hb->r_critical * 1.5};
        int lc_rc = hopf_detect_lc(laser_ode, x_lc, p_lc, 2, hb->r_critical * 0.5,
                                     0.01, 80, 300, &amp, &period);
        if (lc_rc >= 0) {
            printf("  Limit cycle: amp=%.6f period=%.6f\n", amp, period);
        }

        /* Print intensity vs pump curve */
        printf("  Intensity vs pump:\n");
        for (int i = 0; i <= 8; i++) {
            double Pp = P_range[0] + (P_range[1] - P_range[0]) * i / 8.0;
            double I = (Pp > *P_th) ? (Pp - *P_th) / (kappa * gamma / g) : 0.0;
            printf("    P=%.3f: I=%.6f %s\n", Pp, I, Pp > *P_th ? "(lasing)" : "(off)");
        }
    }
    hopf_free(hb);
    return rc;
}

int main(void) {
    double kappa = 10.0, gamma = 1.0, g = 0.5;
    double P_range[] = {0.5, 5.0};
    double P_th, l1;
    int rc = bif_laser_threshold(kappa, gamma, g, P_range, &P_th, &l1);
    printf("=== bif_app2: Laser Hopf Threshold ===\n");
    printf("kappa=%.1f gamma=%.1f g=%.1f -> P_th=%.4f\n", kappa, gamma, g, P_th);
    return rc >= 0 ? 0 : 1;
}
