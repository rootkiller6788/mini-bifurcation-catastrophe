/* ncont_app2.c -- L7-2: Chemical reactor multiplicity (Uppal 1974).
 * CSTR with exothermic reaction: dT/dt = (T0-T)/tau + B*exp(-Ea/RT)*C.
 * Exhibits multiple steady states (ignition/extinction hysteresis loop)
 * as the Damkohler number (Da) or feed temperature varies. */
#include "../include/ncont_arclength.h"
#include "../include/ncont_bifurcation.h"
#include <math.h>
#include <stdio.h>

static void cstr_sys(const double* x, int n, double lambda, void* p, double* r) {
    double* pars = (double*)p;
    double T0 = pars[0], tau = pars[1], B = pars[2], EaR = pars[3];
    double T = x[0] + 300.0;  /* scale for numerics */
    r[0] = (T0 - T) / tau + B * exp(-EaR / (T + 273.15)) * lambda;
}
static void cstr_jac(const double* x, int n, double lambda, void* p, double* J) {
    double* pars = (double*)p;
    double tau = pars[1], B = pars[2], EaR = pars[3];
    double T = x[0] + 300.0;
    J[0] = -1.0 / tau + B * lambda * exp(-EaR/(T+273.15)) * EaR/((T+273.15)*(T+273.15));
}

int ncont_reactor_hysteresis(double Da_range, double* ignition, double* extinction) {
    double pars[] = {300.0, 10.0, 5000.0, 8000.0};  /* T0, tau, B, Ea/R */
    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.05; cfg.ds_max = 0.2; cfg.ds_min = 0.001;
    cfg.max_points = 120; cfg.detect_bifurcations = true;
    cfg.adaptive_stepsize = true;

    double x0[] = {20.0};  /* T - 300, starting on lower branch */
    NCONT_State* state = NULL;

    NCONT_Result res = ncont_continue(cstr_sys, cstr_jac, x0, 1, 0.05, &cfg, pars, &state);

    printf("[ncont_app2] CSTR multiplicity: Da_max=%.2f\n", Da_range);
    printf("  Points: %d, bif: %d, exit: %d\n",
           res.total_points, res.n_bifurcations_detected, res.exit_flag);

    /* Find ignition and extinction points (limit points) */
    NCONT_Bifurcation bif_list[5];
    int nb = ncont_scan_bifurcations(state, cstr_jac, 1, pars, bif_list, 5);

    int ign_idx = -1, ext_idx = -1;
    for (int i = 0; i < nb; i++) {
        if (bif_list[i].type == NCONT_BIF_LIMIT_POINT) {
            if (ign_idx < 0) { ign_idx = i; *ignition = bif_list[i].lambda; }
            else { ext_idx = i; *extinction = bif_list[i].lambda; }
        }
    }

    printf("  Ignition point (Da): %.6f\n", *ignition);
    printf("  Extinction point (Da): %.6f\n", *extinction);
    printf("  Hysteresis window: Da=[%.6f, %.6f]\n", *extinction, *ignition);

    /* Print branch */
    for (int i = 0; i < state->n_points; i += state->n_points/6) {
        printf("  Da=%.4f T=%.1f\n", state->branch[i].lambda,
               state->branch[i].x[0] + 300.0);
    }

    ncont_state_free(state);
    return (ign_idx >= 0 && ext_idx >= 0) ? 0 : 1;
}

int main(void) {
    double ig, ex;
    ncont_reactor_hysteresis(1.0, &ig, &ex);
    printf("=== ncont_app2: CSTR Hysteresis ===\n");
    printf("ignition=%.4f extinction=%.4f\n", ig, ex);
    return 0;
}
