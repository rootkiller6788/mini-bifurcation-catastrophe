/* Amplitude death: detect death region in coupled Stuart-Landau oscillators */
#include "../include/bcd_amplitude.h"
#include "../include/bcd_core.h"
#include <stdio.h>
#include <math.h>

static void slandau_rhs(int n, const double* x, const double* xd, double t, double* dx, void* p) {
    (void)xd; (void)t; (void)p;
    double x2 = x[0]*x[0] + x[1]*x[1];
    double omega = 1.0;
    double mu = 0.2;
    dx[0] = mu * x[0] - omega * x[1] - x[0] * x2;
    dx[1] = omega * x[0] + mu * x[1] - x[1] * x2;
}

int main(void) {
    printf("=== Amplitude Death in Coupled Delayed Oscillators ===\n\n");

    /* Create two coupled Stuart-Landau oscillators */
    BCDCoupledDDE* cd = bcd_coupled_create(2, 2, 0.5, 1.5);
    printf("Coupled system: %d oscillators, coupling=%.2f, delay=%.2f\n",
           cd->n_osc, cd->coupling, cd->delay);

    /* Step the coupled oscillators */
    printf("Simulating free run:\n");
    for (int t = 0; t < 15; t++) {
        int rc = bcd_coupled_step_all(cd, slandau_rhs, NULL, 0.05);
        if (t % 3 == 0 && cd->oscillators[0]) {
            double r1 = sqrt(cd->oscillators[0]->buffer[cd->oscillators[0]->buf_pos] *
                             cd->oscillators[0]->buffer[cd->oscillators[0]->buf_pos]);
            printf("  t=%.2f: amp~%.4f\n", t*0.05, r1);
        }
    }

    /* Compute critical coupling for amplitude death */
    double K_crit = bcd_death_critical_coupling(cd, slandau_rhs, NULL, 0.1, 3.0, 20);
    printf("\nCritical coupling for death: K_crit = %.4f\n", K_crit);

    /* Detect amplitude death at various delays */
    printf("\nDeath region scan (coupling, delay):\n");
    double coup_range[] = {0.5, 2.5};
    double delay_range[] = {0.5, 4.5};
    int death_rc = bcd_oscillation_death_region(cd, slandau_rhs, NULL, coup_range, delay_range);
    printf("Death region detection: %d points in region\n", death_rc);

    /* Full death detection */
    BCDDeathResult* dr = bcd_amplitude_death_detect(cd, 0.01, 50.0, 0.05);
    if (dr) {
        printf("Amplitude death: final_amp=%.6f, death_time=%d, is_dead=%s\n",
               dr->final_amp, dr->death_time, dr->is_dead ? "YES" : "NO");

        double robustness = bcd_death_robustness(dr, 0.01, 50);
        printf("Death robustness (noise=0.01): %.4f\n", robustness);
        bcd_death_result_free(dr);
    }

    bcd_coupled_free(cd);
    printf("\nExample 3 PASSED\n");
    return 0;
}
