#include "bifurcation_pde.h"
#include "numerical_pde.h"
#include "reaction_diffusion.h"
#include <stdio.h>
#include <math.h>

/* Example 2: Bifurcation analysis in PDEs.
 * Track the primary bifurcation from homogeneous state
 * to patterned state as a parameter is varied. */
int main(void) {
    printf("=== PDE Pattern Example 2: Bifurcation Scan ===
");
    int N = 128;
    double dx = 0.2, Du = 1.0, Dv = 25.0;
    RDField* f = rd_create(N, dx, Du, Dv);

    /* Scan over activator production rate alpha */
    printf("Scanning alpha parameter:
");
    printf("  alpha    amplitude    wavelength
");
    for (int ia = 0; ia <= 10; ia++) {
        double alpha = 0.1 + 0.1 * ia;
        TuringParams tp = {alpha, 1.0, -1.0, 0.0, Du, Dv};
        /* Reset field */
        for (int i = 0; i < N; i++) {
            f->u[i] = 1.0 + 0.005*((double)rand()/RAND_MAX-0.5);
            f->v[i] = 0.5;
        }
        /* Simulate */
        for (int s = 0; s < 2000; s++)
            rd_step_euler_1d(f, 0.001, tp);
        double amp = rd_pattern_amplitude(f);
        double wl = rd_pattern_wavelength(f);
        printf("  %.2f    %.6f    %.2f
", alpha, amp, wl);
    }
    rd_free(f);
    return 0;
}
