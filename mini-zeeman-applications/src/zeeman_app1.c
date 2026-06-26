/* zeeman_app1.c -- L7: Heartbeat arrhythmia cusp (Zeeman 1977).
 * Zeeman's heartbeat model: cusp catastrophe on the (tension, fiber_length)
 * plane. The heartbeat cycle traverses the cusp region, with systole and
 * diastole as the two stable branches. Arrhythmias correspond to
 * pathological changes in the cusp control parameters. */
#include "../include/zeeman_core.h"
#include "../include/zeeman_biology.h"
#include "../include/zeeman_dynamics.h"
#include <math.h>
#include <stdio.h>

int zeeman_heart_analysis(double systolic, double diastolic,
                          double* beat_period, int* arrhythmia_risk,
                          double* ejection_fraction) {
    /* Initialize Zeeman cusp model */
    CuspCatastrophe* cusp = cc_create();
    double a = -1.0;  /* normal factor */
    double b = 0.0;   /* splitting factor */

    printf("[zeeman_app1] Heartbeat analysis: systolic=%.1f diastolic=%.1f\n",
           systolic, diastolic);

    /* Map blood pressure to cusp parameters */
    double pressure_diff = systolic - diastolic;
    a = -1.0 + 0.5 * (pressure_diff - 40.0) / 40.0;
    b = 0.2 * (diastolic - 80.0) / 20.0;

    /* Cusp bifurcation set check */
    cusp_bifurcation_set_point(a, &a, &b);
    bool near_bif = fabs(4.0*a*a*a + 27.0*b*b) < 0.01;

    /* Compute equilibrium manifold */
    double eq_upper[50], eq_lower[50], eq_unst[50];
    int n_upper, n_lower, n_unst;
    zeeman_equilibrium_manifold_section(-2.0, 2.0, b, 40,
        eq_upper, eq_lower, eq_unst, &n_upper, &n_lower, &n_unst);

    *beat_period = 0.8 + 0.1 * fabs(a);
    *arrhythmia_risk = near_bif ? 2 : (fabs(a) > 1.5 ? 1 : 0);
    *ejection_fraction = 65.0 - 5.0 * fabs(a);

    printf("  Cusp params: a=%.3f b=%.3f\n", a, b);
    printf("  Near bifurcation set: %s\n", near_bif ? "YES (critical!)" : "no");
    printf("  Equilibrium sections: upper=%d lower=%d unst=%d\n",
           n_upper, n_lower, n_unst);
    printf("  Beat period: %.3f s\n", *beat_period);
    printf("  Arrhythmia risk: %s\n",
           *arrhythmia_risk >= 2 ? "HIGH" : *arrhythmia_risk >= 1 ? "MODERATE" : "LOW");
    printf("  Ejection fraction: %.1f%%\n", *ejection_fraction);

    /* Heartbeat model dynamics */
    HeartbeatModel* hb = (HeartbeatModel*)calloc(1, sizeof(HeartbeatModel));
    hb->fiber_length = 1.0;
    hb->tension = 0.5;
    printf("  Cardiac cycle (10 steps):\n");
    for (int i = 0; i < 10; i++) {
        double dx, dy;
        heartbeat_rhs(hb, &dx, &dy);
        hb->fiber_length += 0.02 * dx;
        hb->tension += 0.02 * dy;
        printf("    t=%.2f: length=%.4f tension=%.4f\n", i*0.05, hb->fiber_length, hb->tension);
    }
    free(hb);

    /* Detect bifurcation crossings in parameter sweep */
    double a_path[] = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0};
    double b_path[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    bool crossings[10];
    int idx;
    bool has_cross = zeeman_detect_bifurcation_crossing(a_path, b_path, 7, &idx);
    printf("  Bifurcation crossing in parameter path: %s (idx=%d)\n",
           has_cross ? "YES" : "NO", idx);

    cc_free(cusp);
    return 0;
}

int main(void) {
    double period, ef;
    int risk;
    zeeman_heart_analysis(120.0, 80.0, &period, &risk, &ef);
    printf("=== zeeman_app1: Heartbeat Cusp ===\n");
    printf("period=%.3f risk=%d ef=%.1f%%\n", period, risk, ef);
    return 0;
}
