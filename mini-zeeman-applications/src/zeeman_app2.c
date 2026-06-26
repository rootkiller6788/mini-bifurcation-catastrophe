/* zeeman_app2.c -- L7-2: Stock market crash cusp (Zeeman 1974).
 * Zeeman's cusp catastrophe model of financial markets:
 * - Normal factor (a): market sentiment / fundamental value
 * - Splitting factor (b): excess demand / speculation
 * The cusp region produces bimodal behavior: bull and bear markets
 * coexist. A crash is a catastrophic jump from the upper sheet
 * (bull) to the lower sheet (bear) at the cusp boundary. */
#include "../include/zeeman_core.h"
#include "../include/zeeman_dynamics.h"
#include <math.h>
#include <stdio.h>

int zeeman_market_crash(double price_index, double speculation_level,
                        double volatility, double* crash_probability,
                        int* crash_risk_level, double* recovery_time) {
    /* Map financial variables to cusp control parameters */
    double a = (price_index - 1000.0) / 500.0;  /* normal factor: deviation from "fair" */
    double b = speculation_level - 0.5;         /* splitting factor: excess speculation */

    printf("[zeeman_app2] Market crash analysis:\n");
    printf("  Price: %.1f, speculation: %.2f, volatility: %.2f\n",
           price_index, speculation_level, volatility);

    /* Cusp bifurcation set: 4a^3 + 27b^2 = 0 */
    double bif_dist = 4.0 * a * a * a + 27.0 * b * b;
    double dist_to_cusp = fabs(bif_dist);

    printf("  Cusp params (a,b)=(%.4f,%.4f) dist_to_bif_set=%.6f\n", a, b, dist_to_cusp);

    /* Crash probability: increases near bifurcation set and with volatility */
    *crash_probability = volatility / (dist_to_cusp + 0.001);
    if (*crash_probability > 1.0) *crash_probability = 1.0;

    /* Risk level classification */
    if (dist_to_cusp < 0.01 || *crash_probability > 0.7)
        *crash_risk_level = 3;  /* CRITICAL: near cusp boundary */
    else if (*crash_probability > 0.3)
        *crash_risk_level = 2;  /* HIGH: approaching cusp */
    else if (b > 0.2 || a < -1.0)
        *crash_risk_level = 1;  /* MODERATE: in cusp region */
    else
        *crash_risk_level = 0;  /* LOW: stable */

    /* Recovery time: depends on how far into bear territory */
    *recovery_time = (a < -1.5) ? 24.0 : (a < -0.5) ? 12.0 : 6.0;

    const char* levels[] = {"LOW", "MODERATE", "HIGH", "CRITICAL"};
    printf("  Crash probability: %.2f%%\n", (*crash_probability) * 100.0);
    printf("  Risk level: %s\n", levels[*crash_risk_level]);
    printf("  Estimated recovery: %.1f months\n", *recovery_time);

    /* Equilibrium manifold analysis */
    CuspCatastrophe* cusp = cc_create();
    cusp->a = &a; cusp->b = &b;

    double eq_up[50], eq_lo[50], eq_un[50];
    int nu, nl, nus;
    zeeman_equilibrium_manifold_section(-3.0, 1.0, b, 30,
        eq_up, eq_lo, eq_un, &nu, &nl, &nus);

    printf("  Market equilibrium sheets: bull=%d pts bear=%d pts unst=%d pts\n", nu, nl, nus);

    /* Hysteresis: once in bear market, harder to recover */
    bool in_bear = a < 0 && b > 0;
    if (in_bear) {
        double new_a = a + 0.5;
        double new_dist = 4.0*new_a*new_a*new_a + 27.0*b*b;
        printf("  Hysteresis: need a=%.2f to exit bear (current a=%.2f, dist=%.4f)\n",
               new_a, a, new_dist);
    }

    /* Detect bifurcation crossing in hypothetical market trajectory */
    double a_traj[] = {-2.0, -1.5, -1.0, -1.2, -0.8};
    double b_traj[] = {0.3, 0.3, 0.4, 0.4, 0.5};
    int cross_idx;
    bool crossed = zeeman_detect_bifurcation_crossing(a_traj, b_traj, 5, &cross_idx);
    printf("  Market trajectory crosses bifurcation: %s\n", crossed ? "YES (CRASH!)" : "no");

    cc_free(cusp);
    return 0;
}

int main(void) {
    double prob, recovery;
    int risk;
    zeeman_market_crash(850.0, 0.7, 0.3, &prob, &risk, &recovery);
    printf("=== zeeman_app2: Stock Market Cusp ===\n");
    printf("crash_prob=%.2f risk=%d recovery=%.1f months\n", prob, risk, recovery);
    return 0;
}
