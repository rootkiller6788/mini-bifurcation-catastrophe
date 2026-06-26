#include "reaction_diffusion.h"
#include "pattern_formation.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Example 1: Turing pattern formation in 1D.
 * Demonstrates how diffusion-driven instability creates
 * spontaneous spatial patterns from a homogeneous state. */
int main(void) {
    printf("=== PDE Pattern Example 1: Turing Instability ===
");
    int N = 256;
    double dx = 0.1, Du = 1.0, Dv = 40.0;
    RDField* f = rd_create(N, dx, Du, Dv);
    if (!f) { printf("Failed to create field
"); return 1; }

    /* Set initial: small random perturbation around equilibrium */
    TuringParams tp = {0.5, 1.0, -1.0, 0.0, Du, Dv};
    for (int i = 0; i < N; i++) {
        f->u[i] = 1.0 + 0.01 * ((double)rand()/RAND_MAX - 0.5);
        f->v[i] = 0.5 + 0.01 * ((double)rand()/RAND_MAX - 0.5);
    }

    /* Check Turing conditions */
    int can_form;
    double k2_min, k2_max;
    rd_turing_conditions(tp, &k2_min, &k2_max, &can_form);
    printf("Turing possible: %s (k2 in [%.3f, %.3f])
",
           can_form ? "YES" : "no", k2_min, k2_max);

    /* Simulate a few steps */
    double dt = 0.001;
    for (int step = 0; step < 1000; step++)
        rd_step_euler_1d(f, dt, tp);

    /* Measure pattern */
    double amp = rd_pattern_amplitude(f);
    double wl = rd_pattern_wavelength(f);
    printf("Pattern amplitude: %.6f
", amp);
    printf("Pattern wavelength: %.3f
", wl);

    rd_free(f);
    return 0;
}
