/* Zeeman applications: heartbeat, buckling, and regime shift */
#include "../include/tec_applications.h"
#include "../include/tec_bifurcation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Zeeman Catastrophe Applications ===\n\n");

    /* --- Heartbeat model --- */
    printf("Heartbeat model (Zeeman 1977):\n");
    TEC_HeartbeatState hb = tec_heartbeat_init();
    printf("  Init: tension=%.4f, fiber_length=%.4f\n", hb.tension, hb.fiber_length);
    for (int i = 0; i < 10; i++) {
        tec_heartbeat_step(&hb, 0.05);
        double ecg = tec_heartbeat_get_ecg(&hb);
        printf("  t=%.2f: ecg=%.4f tension=%.4f length=%.4f\n",
               (i+1)*0.05, ecg, hb.tension, hb.fiber_length);
    }

    /* --- Euler buckling --- */
    printf("\nEuler Buckling (fold catastrophe):\n");
    double E = 200e9, I = 1e-6, L = 2.0;
    double P_crit = tec_buckling_critical_load(E, I, L);
    printf("  E=%.1e I=%.1e L=%.1f => P_crit = %.2f N\n", E, I, L, P_crit);
    for (double ratio = 0.5; ratio <= 1.5; ratio += 0.25) {
        double P = ratio * P_crit;
        double defl = tec_buckling_deflection(P, P_crit, 0.01);
        printf("  P/P_crit=%.2f deflection=%.6f m\n", ratio, defl);
    }

    /* --- Ecological regime shift --- */
    printf("\nEcological regime shift (fold-cusp):\n");
    double K = 100.0, h = 15.0, r = 0.5;
    double threshold = tec_regime_shift_threshold(K, h, r);
    printf("  K=%.1f, h=%.1f, r=%.2f => threshold=%.4f\n", K, h, r, threshold);

    /* --- Ship stability (butterfly) --- */
    printf("\nShip stability (butterfly catastrophe):\n");
    double bfly[] = {-2.0, -0.5, 0.3, 0.0};
    double capsize = tec_ship_capsize_angle(bfly);
    printf("  Capsize angle: %.2f degrees\n", capsize);
    for (double deg = 0; deg <= 40; deg += 10) {
        double rad = deg * 3.141592653589793 / 180.0;
        double moment = tec_ship_restoring_moment(rad, bfly[0], bfly[1], bfly[2], bfly[3]);
        printf("  heel=%.0f deg: restoring moment=%.4f\n", deg, moment);
    }

    /* --- Phase transition (Landau) --- */
    printf("\nPhase transition (Landau theory):\n");
    double Tc = 100.0, a4 = -1.0, a6 = 0.1;
    for (double T = 80; T <= 120; T += 10) {
        double eq = tec_landau_equilibrium(T, Tc, a4, a6);
        double F = tec_landau_free_energy(eq, T, Tc, a4, a6);
        printf("  T=%.0f: eq=%.4f, F=%.4f\n", T, eq, F);
    }

    printf("\nExample 3 PASSED\n");
    return 0;
}
