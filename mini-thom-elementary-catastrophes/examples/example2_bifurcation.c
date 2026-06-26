/* Thom catastrophes: bifurcation diagrams and hysteresis analysis */
#include "../include/tec_bifurcation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Bifurcation Diagrams & Hysteresis ===\n\n");

    /* --- Fold bifurcation diagram --- */
    TEC_BifurcationDiagram* fold_bd = tec_bif_diagram_fold(30, 2.0);
    printf("Fold diagram: %d parameter values\n", fold_bd->n_params);
    printf("Equilibrium manifold (x vs a):\n");
    for (int i = 0; i < fold_bd->n_params; i++) {
        printf("  a=%.2f: %d equilibria", fold_bd->param_range[i], fold_bd->n_eq_per_param[i]);
        for (int j = 0; j < fold_bd->n_eq_per_param[i]; j++) {
            printf(" x=%.3f", fold_bd->equilibria[i * fold_bd->max_eq + j]);
        }
        printf("\n");
    }
    tec_bif_diagram_free(fold_bd);

    /* --- Cusp bifurcation diagram (fixed a, vary b) --- */
    printf("\n--- Cusp slice: a=-1.5 (inside cusp) ---\n");
    TEC_BifurcationDiagram* cusp_bd = tec_bif_diagram_cusp(20, 1.5, -1.5);
    for (int i = 0; i < cusp_bd->n_params; i++) {
        printf("  b=%.2f: %d equilibria", cusp_bd->param_range[i], cusp_bd->n_eq_per_param[i]);
        for (int j = 0; j < cusp_bd->n_eq_per_param[i]; j++) {
            printf(" x=%.3f", cusp_bd->equilibria[i * cusp_bd->max_eq + j]);
        }
        printf("\n");
    }
    tec_bif_diagram_free(cusp_bd);

    /* --- Hysteresis regions --- */
    printf("\n--- Hysteresis regions in cusp ---\n");
    double test_pts[][2] = {{-2, 0}, {-1, 0.5}, {-0.5, 0}, {0.5, 0}};
    for (int i = 0; i < 4; i++) {
        bool h = tec_has_hysteresis(TEC_CUSP, test_pts[i]);
        printf("  (a,b)=(%.1f,%.1f): hysteresis=%s", test_pts[i][0], test_pts[i][1], h ? "YES" : "NO");
        if (h) printf(", width=%.4f", tec_hysteresis_loop_width(TEC_CUSP, test_pts[i]));
        printf("\n");
    }

    /* --- Delay convention jump --- */
    printf("\n--- Delay convention ---\n");
    double jump = tec_delay_convention_jump(TEC_CUSP, test_pts[0], -2.0, -0.1);
    printf("Jump at a=%.1f,b=%.1f: jump=%.4f\n", test_pts[0][0], test_pts[0][1], jump);

    /* --- Maxwell set for cusp --- */
    TEC_MaxwellPoint* mw = tec_maxwell_set_cusp(10, 3.0);
    if (mw) {
        printf("\nMaxwell set (equal-depth minima):\n");
        for (int i = 0; i < 3 && i < 10; i++) {
            printf("  a=%.3f: x1=%.4f x2=%.4f V_min=%.4f\n",
                   mw[i].control[0], mw[i].x1, mw[i].x2, mw[i].V_min);
        }
        tec_maxwell_free(mw, 10);
    }

    printf("\nExample 2 PASSED\n");
    return 0;
}
