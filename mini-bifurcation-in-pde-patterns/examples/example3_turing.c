#include "reaction_diffusion.h"
#include "pattern_formation.h"
#include <stdio.h>
#include <math.h>

/* Example 3: Dispersion relation and mode selection.
 * Compute the growth rate Re(lambda(k^2)) for each
 * wavenumber k to predict the dominant pattern mode. */
int main(void) {
    printf("=== PDE Pattern Example 3: Dispersion Relation ===
");
    TuringParams tp = {0.5, 1.0, -1.0, 0.0, 1.0, 40.0};

    /* Compute dispersion for range of k^2 */
    printf("  k^2       Re(lambda)
");
    double max_growth = -1e100;
    double best_k2 = 0.0;
    for (int i = 0; i <= 100; i++) {
        double k2 = 0.01 * (double)i;
        double growth = rd_dispersion_relation(k2, tp);
        printf("  %.4f    %+.6f
", k2, growth);
        if (growth > max_growth) {
            max_growth = growth;
            best_k2 = k2;
        }
    }
    printf("Most unstable mode: k^2 = %.4f (wavelength = %.2f)
",
           best_k2, 2.0*M_PI/sqrt(best_k2));

    /* Check Turing conditions */
    int can_form;
    double k2_min, k2_max;
    rd_turing_conditions(tp, &k2_min, &k2_max, &can_form);
    printf("Turing band: k^2 in [%.4f, %.4f]
", k2_min, k2_max);

    return 0;
}
