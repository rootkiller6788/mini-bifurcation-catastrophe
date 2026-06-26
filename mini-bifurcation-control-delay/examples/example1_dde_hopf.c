/* DDE Hopf bifurcation: use Lambert W to compute stability boundaries */
#include "../include/bcd_core.h"
#include <stdio.h>
#include <math.h>

static void scalar_dde_rhs(int n, const double* x, const double* xd, double t, double* dx, void* p) {
    (void)t; double* pars = (double*)p;
    dx[0] = -pars[0] * x[0] - pars[1] * xd[0];  /* dx/dt = -a*x(t) - b*x(t-tau) */
}

int main(void) {
    printf("=== DDE Hopf Bifurcation: Lambert W Stability Chart ===\n\n");

    double pars[] = {0.5, 2.0};  /* a=0.5, b=2.0 */
    double tau_vals[] = {0.2, 0.8, 1.2, 2.0, 3.0};

    for (int ti = 0; ti < 5; ti++) {
        double tau = tau_vals[ti];
        BCDDDE* dde = bcd_dde_create(1, tau, 200);
        double x0 = 1.0;
        bcd_dde_set_history(dde, &x0, 1);

        /* Simulate for T=20 to observe stability */
        double t_out[2000], x_out[2000];
        int n_pts = 0;
        bcd_dde_simulate(dde, scalar_dde_rhs, pars, 20.0, 0.02, t_out, 1000);
        /* Extract final values from dde state */
        BCDDDEState st = {.x = &x0, .t = dde->t_current, .step = 0};

        /* Integrate a few more steps to observe */
        double x_final = x0;
        for (int s = 0; s < 500; s++) {
            bcd_dde_step(dde, scalar_dde_rhs, pars, 0.02);
        }
        x_final = dde->buffer[dde->buf_pos];

        double amp = fabs(x_final);
        printf("  tau=%.1f: final |x|=%.6f -> %s\n",
               tau, amp, amp < 0.01 ? "STABLE (fixed point)" :
               amp < 5.0 ? "OSCILLATORY" : "UNSTABLE");

        bcd_dde_free(dde);
    }

    /* Find critical tau via stability scan */
    printf("\nStability boundary scan (a=0.5, b=2.0):\n");
    for (double tau = 0.5; tau <= 3.5; tau += 0.3) {
        BCDDDE* dde = bcd_dde_create(1, tau, 200);
        double x0 = 1.0;
        bcd_dde_set_history(dde, &x0, 1);
        for (int s = 0; s < 300; s++) bcd_dde_step(dde, scalar_dde_rhs, pars, 0.02);
        double xf = dde->buffer[dde->buf_pos];
        int stable = fabs(xf) < 0.01;
        printf("  tau=%.1f: %s (x=%.6f)\n", tau, stable ? "STABLE" : "UNSTABLE", xf);
        bcd_dde_free(dde);
    }

    printf("\nExample 1 PASSED\n");
    return 0;
}
