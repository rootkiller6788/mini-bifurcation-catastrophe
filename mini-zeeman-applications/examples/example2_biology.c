#include "zeeman_biology.h"
#include "zeeman_dynamics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=======================================\n");
    printf("  Zeeman Biological Applications — Demo 2\n");
    printf("  Heartbeat, Nerve Impulse, Cell Fate\n");
    printf("=======================================\n\n");

    /* 1. Heartbeat simulation */
    printf("--- 1. Heartbeat Model (Zeeman 1972) ---\n");
    HeartbeatModel* hb = heartbeat_create(-3.0, 0.0, 0.08, 0.3);
    printf("Initial: "); heartbeat_print(hb);

    printf("\nSimulating 30 time units...\n");
    heartbeat_simulate(hb, 30.0, 0.01);
    printf("Final: "); heartbeat_print(hb);

    double period = heartbeat_period(hb, 20.0, 0.01);
    printf("Estimated period: %.3f\n", period);

    heartbeat_free(hb);

    /* 2. Nerve impulse */
    printf("\n--- 2. Nerve Impulse Model (Zeeman 1973) ---\n");
    /* Subthreshold stimulus */
    NerveImpulse* n_sub = nerve_impulse_create(-1.0, 0.5, 0.1);
    printf("Subthreshold (I=0.1):\n");
    nerve_impulse_simulate(n_sub, 10.0, 0.01);
    int spikes_sub = nerve_impulse_spike_count(n_sub, 10.0, 0.01);
    printf("  Spikes: %d, Firing: %s\n",
           spikes_sub, nerve_impulse_is_firing(n_sub) ? "YES" : "NO");
    nerve_impulse_free(n_sub);

    /* Superthreshold stimulus */
    NerveImpulse* n_super = nerve_impulse_create(-1.0, 0.5, 1.0);
    printf("\nSuperthreshold (I=1.0):\n");
    printf("  Threshold: %.4f\n", nerve_impulse_threshold(n_super));
    nerve_impulse_simulate(n_super, 10.0, 0.01);
    int spikes_super = nerve_impulse_spike_count(n_super, 10.0, 0.01);
    printf("  Spikes: %d, Firing: %s\n",
           spikes_super, nerve_impulse_is_firing(n_super) ? "YES" : "NO");
    nerve_impulse_free(n_super);

    /* 3. Cell differentiation */
    printf("\n--- 3. Cell Differentiation (Thom 1975) ---\n");
    printf("Varying morphogen from 0 to 2.0:\n");

    for (double m = 0.0; m <= 2.0; m += 0.4) {
        CellDiffModel* cd = cell_diff_create(m, 0.5, 2.0);
        cell_diff_simulate(cd, 5.0, 0.01);
        const char* fate = (cd->x > 0.5) ? "DIFFERENTIATED"
                         : (cd->x < -0.5) ? "ALTERNATIVE"
                         : "UNDIFFERENTIATED";
        printf("  morphogen=%.1f  x=%.4f  stable_states=%d  fate=%s\n",
               m, cd->x, cell_diff_stable_states(cd), fate);
        cell_diff_free(cd);
    }

    /* 4. Slow-fast relaxation oscillations */
    printf("\n--- 4. Slow-Fast Cusp Relaxation Oscillations ---\n");
    SlowFastCusp* sf = slowfast_cusp_create(0.05, -2.0, 0.0, 3.0, 0.5, 1.0);
    printf("Initial: "); slowfast_cusp_print(sf);

    printf("Simulating 30 time units...\n");
    slowfast_cusp_simulate(sf, 30.0, 0.01);
    printf("Final: "); slowfast_cusp_print(sf);

    double sf_period = slowfast_cusp_period(sf, 30.0, 0.01);
    printf("Relaxation oscillation period: %.3f\n", sf_period);
    slowfast_cusp_free(sf);

    printf("\n=======================================\n");
    printf("  Demo complete.\n");
    printf("=======================================\n");
    return 0;
}