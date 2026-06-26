/* Hopf bifurcation: detect and classify in a 2D system */
#include "../include/hopf_bifurcation.h"
#include <stdio.h>
#include <math.h>

static void hopf_test_sys(const double* x, double* params, double* dx, int n) {
    (void)n;
    double mu = params[0], omega = params[1];
    double x2 = x[0]*x[0] + x[1]*x[1];
    dx[0] = mu * x[0] - omega * x[1] - x[0] * x2;  /* supercritical Hopf normal form */
    dx[1] = omega * x[0] + mu * x[1] - x[1] * x2;
}

int main(void) {
    printf("=== Hopf Bifurcation Detection ===\n\n");

    double x0[] = {0.1, 0.1};
    double params[] = {-0.2, 1.5};  /* mu=-0.2, stable focus */
    HopfBifurcation* hb = hopf_create();

    int rc = hopf_detect(hopf_test_sys, x0, params, 2, 0, -1.0, 1.0, 0.05, 50, 1e-8, hb);
    if (rc < 0) { printf("Detection FAILED\n"); hopf_free(hb); return 1; }

    printf("Detected at mu_critical = %.6f\n", hb->r_critical);
    printf("Critical frequency omega = %.6f\n", hb->omega_critical);
    printf("First Lyapunov coefficient l1 = %.6f\n", hb->first_lyapunov);
    printf("Type: %s\n", hb->is_supercritical ? "SUPERCRITICAL (stable LC born)" : "SUBCRITICAL (unstable LC)");

    /* Compute Lyapunov coefficient directly */
    double l1 = hopf_first_lyapunov(hopf_test_sys, x0, params, 2, 0, hb->omega_critical, 1e-6);
    printf("Computed l1 %.6f => %s\n", l1, l1 < 0 ? "Supercritical" : "Subcritical");

    /* Detect limit cycle amplitude */
    double amp, period;
    double x_lc[] = {0.2, 0.0};
    double p_lc[] = {0.3, 1.5};
    int lc_rc = hopf_detect_lc(hopf_test_sys, x_lc, p_lc, 2, 0.3, 0.01, 100, 500, &amp, &period);
    if (lc_rc >= 0) {
        printf("Limit cycle: amplitude=%.6f, period=%.6f\n", amp, period);
    }

    hopf_print(hb);
    hopf_free(hb);
    printf("\nExample 2 PASSED\n");
    return 0;
}
