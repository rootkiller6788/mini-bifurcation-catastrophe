/* tec_app2.c -- L7-2: Phase transition cusp (Landau theory, Van der Waals).
 * The Van der Waals equation of state P = RT/(V-b) - a/V^2 exhibits a
 * first-order liquid-gas phase transition below the critical temperature.
 * Near the critical point (Tc, Pc, Vc), the free energy is a cusp catastrophe
 * with order parameter (V-Vc) and control parameters (T-Tc, P-Pc). */
#include "../include/tec_core.h"
#include "../include/tec_bifurcation.h"
#include <math.h>
#include <stdio.h>

/* Van der Waals parameters for CO2 */
#define VDW_a 0.364  /* J*m^3/mol^2 */
#define VDW_b 4.27e-5 /* m^3/mol */
#define VDW_R 8.314   /* J/(mol*K) */

static double vdw_pressure(double T, double V) {
    if (V <= VDW_b) return 1e10;
    return VDW_R * T / (V - VDW_b) - VDW_a / (V * V);
}

/* Map Van der Waals to cusp catastrophe near critical point */
int tec_phase_transition(double T, double* V_vals, double* P_vals, int n,
                          int* has_hysteresis) {
    double Tc = 8.0 * VDW_a / (27.0 * VDW_R * VDW_b);
    double Vc = 3.0 * VDW_b;
    double Pc = VDW_a / (27.0 * VDW_b * VDW_b);

    printf("[tec_app2] Van der Waals CO2: Tc=%.2f K, Pc=%.4f Pa, Vc=%.2e m^3\n", Tc, Pc, Vc);

    /* Map to cusp control parameters */
    double a_ctrl = (T - Tc) / Tc;  /* normal factor */
    double b_ctrl = 0.0;

    /* Compute isotherm */
    double V_min = 0.8 * Vc, V_max = 3.0 * Vc;
    double dV = (V_max - V_min) / (double)(n - 1);

    printf("[tec_app2] T=%.1f K (%.2f Tc) isotherm:\n", T, T/Tc);
    for (int i = 0; i < n; i++) {
        V_vals[i] = V_min + i * dV;
        P_vals[i] = vdw_pressure(T, V_vals[i]);
        if (i % (n/8) == 0)
            printf("  V=%.3e P=%.4f\n", V_vals[i], P_vals[i]);
    }

    /* Check hysteresis using catastrophe theory */
    double cusp_params[] = {a_ctrl, b_ctrl};
    *has_hysteresis = tec_has_hysteresis(TEC_CUSP, cusp_params) ? 1 : 0;

    if (*has_hysteresis) {
        double hw = tec_hysteresis_loop_width(TEC_CUSP, cusp_params);
        printf("  Hysteresis: YES (width=%.4f in order param)\n", hw);

        /* Maxwell equal-area construction via Maxwell set */
        TEC_MaxwellPoint* mw = tec_maxwell_set_cusp(10, fabs(a_ctrl));
        if (mw) {
            printf("  Maxwell point: a=%.4f, x1=%.4f x2=%.4f\n",
                   mw[0].control[0], mw[0].x1, mw[0].x2);
            tec_maxwell_free(mw, 10);
        }
    } else {
        printf("  Hysteresis: NO (supercritical, T > Tc)\n");
    }

    /* Landau free energy near critical point */
    double phi = tec_landau_equilibrium(T, Tc, -1.0, 0.05);
    double F = tec_landau_free_energy(phi, T, Tc, -1.0, 0.05);
    printf("  Landau: phi=%.6f F=%.6f\n", phi, F);

    return 0;
}

int main(void) {
    double V[50], P[50];
    int has_hyst;
    tec_phase_transition(280.0, V, P, 30, &has_hyst);  /* below Tc=304K */
    printf("=== tec_app2: Phase Transition Cusp ===\n");
    printf("T/Tc=%.2f, hysteresis=%s\n", 280.0/304.13, has_hyst ? "YES" : "NO");
    return 0;
}
