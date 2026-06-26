/* Melnikov method: compute chaos threshold for forced Duffing */
#include "../include/homoclinic_melnikov.h"
#include "../include/homoclinic_orbit.h"
#include <stdio.h>
#include <math.h>

static void duffing_unperturbed(HOMState x, HOMState* dx, void* p) {
    (void)p;
    dx->x[0] = x.x[1];
    dx->x[1] = x.x[0] - x.x[0] * x.x[0] * x.x[0];  /* xddot - x + x^3 = 0 */
}
static double duffing_forcing(double t, void* p) {
    double* pars = (double*)p;
    return pars[0] * cos(pars[1] * t);  /* gamma*cos(omega*t) */
}

static void generate_homoclinic_orbit(HOMState* orbit, int n) {
    /* Analytical: x(t)=sqrt(2)*sech(t), y(t)=-sqrt(2)*sech(t)*tanh(t) */
    double t_start = -5.0, t_end = 5.0, dt = (t_end - t_start) / (double)(n-1);
    for (int i = 0; i < n; i++) {
        double t = t_start + i * dt;
        double sech_t = 1.0 / cosh(t);
        orbit[i].x[0] = 1.414213562373095 * sech_t;       /* sqrt(2)*sech(t) */
        orbit[i].x[1] = -1.414213562373095 * sech_t * tanh(t);
        orbit[i].dim = 2;
    }
}

int main(void) {
    printf("=== Melnikov Method: Forced Duffing Oscillator ===\n\n");

    double forcing_pars[] = {0.3, 1.2};  /* gamma=0.3, omega=1.2 */
    HOMMelnikov* mel = hom_melnikov_create(duffing_unperturbed, NULL,
                                            duffing_forcing, forcing_pars);

    /* Generate unperturbed homoclinic orbit */
    int n_orbit = 200;
    HOMState* orbit = (HOMState*)calloc((size_t)n_orbit, sizeof(HOMState));
    generate_homoclinic_orbit(orbit, n_orbit);

    /* Compute Melnikov integral */
    double M = hom_melnikov_integral(mel, orbit, n_orbit, 0.0, 0.05);
    printf("Melnikov integral M(t0=0) = %.6f\n", M);

    /* Find zero crossings of Melnikov function */
    int nz = hom_melnikov_find_zeros(mel, orbit, n_orbit, 0.05, -20.0, 20.0, 50);
    printf("Zero crossings found: %d\n", nz);
    if (nz > 0) {
        printf("Zeros at:");
        for (int i = 0; i < nz && i < 8; i++) {
            printf(" t0=%.3f", mel->zero_crossings[i]);
        }
        printf("\n");
    }

    /* Check for transverse intersection */
    bool transverse = hom_melnikov_has_transverse_intersection(mel);
    printf("Transverse intersection: %s\n", transverse ? "YES (chaos!)" : "NO");

    /* Bifurcation threshold */
    double threshold = hom_melnikov_bifurcation_threshold(mel);
    printf("Bifurcation threshold (gamma_crit): %.6f\n", threshold);

    /* Chaos prediction */
    bool chaos = hom_melnikov_chaos_prediction(mel, 0.1);
    printf("Chaos prediction: %s\n", chaos ? "CHAOS expected" : "NO chaos");

    hom_melnikov_print(mel);
    hom_melnikov_free(mel);
    free(orbit);
    printf("\nExample 2 PASSED\n");
    return 0;
}
