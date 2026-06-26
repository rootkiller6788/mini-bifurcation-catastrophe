/* bif_app1.c -- L7: Euler buckling pitchfork (Euler 1744, elastic instability).
 * Models an elastic column under axial compression. The trivial straight
 * configuration loses stability at the critical load via a supercritical
 * pitchfork bifurcation. Post-buckling: two symmetric deflected equilibria. */
#include "../include/pitchfork_bifurcation.h"
#include <math.h>
#include <stdio.h>

/* Euler buckling ODE system: state = (lateral_deflection, rotation)
 * dx/dt = y, dy/dt = P*(pi/L)^2*x/EI - (pi/L)^2*x^3 */
static void euler_ode(const double* x, double* params, double* dx, int n) {
    (void)n;
    double P = params[0];    /* axial load */
    double k = params[1];    /* pi^2*EI/L^2 */
    dx[0] = x[1];
    dx[1] = P * x[0] / k - x[0] * x[0] * x[0];
}

int bif_euler_buckling(double E, double I, double L, double* P_crit, double* post_buckling_curve, int n_pts) {
    double k = M_PI * M_PI * E * I / (L * L);  /* Euler critical load */
    *P_crit = k;

    double x0[] = {0.01, 0.0};
    double params[] = {0.5 * k, k};

    PitchforkBifurcation* pf = pf_create();
    int rc = pf_detect(euler_ode, x0, params, 2, 0, 0.3 * k, 2.0 * k, 0.05 * k, 50, 1e-8, pf);

    if (rc >= 0) {
        printf("[bif_app1] Euler buckling: P_crit=%.2f N (expected %.2f N)\n", pf->r_critical, k);
        printf("  PF coeff = %.6f, supercritical=%s, Z2_sym=%s\n",
               pf->pf_coeff, pf->is_supercritical ? "YES" : "NO",
               pf->has_z2_sym ? "YES" : "NO");

        /* Post-buckling branches */
        int nb = pf_compute_branches(pf, euler_ode, x0, params, 2, 0,
                                      0.5 * k, 2.0 * k, n_pts, 50, 1e-8);
        for (int i = 0; i < nb && i < n_pts; i++) {
            double ratio = pf->param_range[i] / k;
            post_buckling_curve[i] = pf->asym_plus[i];
            printf("  P/P_crit=%.3f: asym+=%.6f asym-=%.6f\n",
                   ratio, pf->asym_plus[i], pf->asym_minus[i]);
        }
    }
    pf_free(pf);
    return rc;
}

int main(void) {
    double E = 200e9, I = 1e-6, L = 2.0;
    double P_crit, curve[20];
    int rc = bif_euler_buckling(E, I, L, &P_crit, curve, 20);
    printf("=== bif_app1: Euler Buckling Pitchfork ===\n");
    printf("E=%.1e Pa, I=%.1e m^4, L=%.1f m\n", E, I, L);
    printf("P_crit=%.2f N, post-buckling computed: %s\n", P_crit, rc >= 0 ? "YES" : "NO");
    return rc >= 0 ? 0 : 1;
}
