/* pde_app2.c -- L7-2: Vegetation pattern (Klausmeier 1999, semi-arid).
 * The Klausmeier model: dW/dt = A - LW - R*W*B^2 + Dw*d2W/dx2
 *                      dB/dt = R*J*W*B^2 - M*B + Db*d2B/dx2
 * where W=water, B=biomass. Turing instability creates banded
 * vegetation patterns (tiger bush) observed in semi-arid ecosystems. */
#include "../include/reaction_diffusion.h"
#include "../include/numerical_pde.h"
#include <math.h>
#include <stdio.h>

int pde_vegetation_pattern(double rainfall_A, double loss_L, double uptake_R,
                            double conversion_J, double mortality_M,
                            double Dw, double Db, double* wavelength,
                            double* pattern_speed, double* biomass_peak) {
    /* Homogeneous steady state */
    double B0 = (rainfall_A * R * conversion_J - M * L) / (L * M);
    double W0 = M / (R * conversion_J * B0);

    if (B0 <= 0) {
        printf("[pde_app2] No vegetation possible (B0=%.4f <= 0)\n", B0);
        return -1;
    }

    /* Jacobian at steady state */
    double J11 = -L - uptake_R * B0 * B0;
    double J12 = -2.0 * uptake_R * W0 * B0;
    double J21 = uptake_R * conversion_J * B0 * B0;
    double J22 = 2.0 * uptake_R * conversion_J * W0 * B0 - mortality_M;

    double det = J11 * J22 - J12 * J21;
    double trace = J11 + J22;

    printf("[pde_app2] Vegetation pattern (Klausmeier): rainfall=%.3f\n", rainfall_A);
    printf("  Steady state: W0=%.4f B0=%.4f\n", W0, B0);
    printf("  Jacobian: trace=%.4f det=%.4f\n", trace, det);
    printf("  Turing condition: Dw/Db=%.4f, trace=%s\n",
           Dw/Db, trace < 0 ? "stable" : "unstable");

    /* Turing instability condition for Dv >> Du */
    double ratio = Db / Dw;
    double turing_threshold = (J22 * J22) / (det - J11 * J22);
    printf("  Required Db/Dw > %.2f (actual %.2f): %s\n",
           turing_threshold, ratio,
           ratio > turing_threshold ? "Turing possible!" : "below threshold");

    /* Dominant wavelength */
    double k2_max = sqrt(det / (Dw * Db));
    *wavelength = 2.0 * M_PI / sqrt(k2_max);
    *pattern_speed = 0.0;  /* stationary for Turing */
    *biomass_peak = B0 * 1.5;

    printf("  Dominant pattern wavelength: %.4f m\n", *wavelength);
    printf("  Pattern type: %s (stationary)\n", *pattern_speed < 0.01 ? "Turing" : "traveling wave");

    /* Parameter sensitivity: find critical rainfall */
    printf("  Critical rainfall scan:\n");
    for (double pa = 0.1; pa <= 0.5; pa += 0.1) {
        double B = (pa * uptake_R * conversion_J - mortality_M * loss_L) / (loss_L * mortality_M);
        double Juu = -loss_L - uptake_R * B * B;
        double Jvv = 2.0 * uptake_R * conversion_J * (M/(R*conversion_J*B)) * B - mortality_M;
        double trace2 = Juu + Jvv;
        printf("    A=%.2f: B0=%.4f trace=%.4f %s\n", pa, B, trace2,
               B > 0 && trace2 < 0 ? "(viable+stable)" : B <= 0 ? "(extinct)" : "(unstable)");
    }

    return 0;
}

int main(void) {
    double wl, speed, peak;
    pde_vegetation_pattern(0.3, 0.1, 1.0, 0.5, 0.15, 0.01, 2.0, &wl, &speed, &peak);
    printf("=== pde_app2: Vegetation Pattern ===\n");
    printf("wavelength=%.4f peak=%.4f\n", wl, peak);
    return 0;
}
