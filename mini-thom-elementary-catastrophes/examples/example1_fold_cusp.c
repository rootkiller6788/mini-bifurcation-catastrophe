/* Thom catastrophes: compute equilibria and bifurcation sets for fold & cusp */
#include "../include/tec_core.h"
#include "../include/tec_bifurcation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Thom Catastrophes: Fold & Cusp Geometry ===\n\n");

    /* --- Fold catastrophe (codim 1) --- */
    TEC_Model* fold = tec_model_create(TEC_FOLD);
    double fold_params[] = {-1.0};
    tec_model_set_control(fold, fold_params, 1);

    int n_eq = tec_find_equilibria(fold);
    printf("Fold V=x^3 + a*x, a=%.1f:\n", fold_params[0]);
    printf("  Equilibria found: %d (expect 2 for a<0, 0 for a>0)\n", n_eq);
    for (int i = 0; i < n_eq; i++) {
        printf("  eq[%d]: x=%.4f, V=%.4f, type=%d\n", i,
               fold->equilibria[i].x, fold->equilibria[i].V,
               fold->equilibria[i].type);
    }

    TEC_BifurcationSet* fold_bs = tec_bif_set_fold(20);
    printf("Fold bifurcation set: %d points (a=0 line)\n", fold_bs->n_points);
    tec_bif_set_free(fold_bs);

    /* --- Cusp catastrophe (codim 2) --- */
    TEC_Model* cusp = tec_model_create(TEC_CUSP);
    double cusp_params[] = {-2.0, 0.0};  /* inside cusp => 3 eq */
    tec_model_set_control(cusp, cusp_params, 2);

    n_eq = tec_find_equilibria(cusp);
    printf("\nCusp V=x^4/4 + a*x^2/2 + b*x, (a,b)=(%.1f,%.1f):\n", cusp_params[0], cusp_params[1]);
    printf("  Equilibria: %d (expect 3 inside cusp)\n", n_eq);
    for (int i = 0; i < n_eq; i++) {
        printf("  eq[%d]: x=%.4f, V=%.4f, hess=%.4f, type=%d\n", i,
               cusp->equilibria[i].x, cusp->equilibria[i].V,
               cusp->equilibria[i].hessian_det, cusp->equilibria[i].type);
    }

    /* Hysteresis check */
    bool hyst = tec_has_hysteresis(TEC_CUSP, cusp_params);
    printf("  Has hysteresis: %s\n", hyst ? "YES" : "NO");

    if (hyst) {
        double hw = tec_hysteresis_loop_width(TEC_CUSP, cusp_params);
        printf("  Hysteresis loop width: %.4f\n", hw);
    }

    /* Cusp bifurcation set */
    TEC_BifurcationSet* cusp_bs = tec_bif_set_cusp(40, 3.0);
    printf("\nCusp bifurcation set: %d points on 4a^3+27b^2=0\n", cusp_bs->n_points);
    for (int i = 0; i < 5 && i < cusp_bs->n_points; i++) {
        printf("  pt[%d]: a=%.4f b=%.4f\n", i,
               cusp_bs->points[i].control[0], cusp_bs->points[i].control[1]);
    }
    tec_bif_set_free(cusp_bs);

    /* Potential evaluation across range */
    printf("\nPotential landscape scan (cusp a=-2,b=0):\n");
    for (double x = -2.5; x <= 2.5; x += 0.5) {
        double V = tec_potential_cusp(x, -2.0, 0.0);
        printf("  x=%+.1f V=%.4f\n", x, V);
    }

    tec_model_free(fold);
    tec_model_free(cusp);
    printf("\nExample 1 PASSED\n");
    return 0;
}
