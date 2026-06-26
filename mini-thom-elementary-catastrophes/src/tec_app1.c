/* tec_app1.c -- L7: Caustics in optics (Berry 1976, fold and cusp catastrophes).
 * In ray optics, caustics are envelopes of light rays where intensity diverges.
 * The fold caustic (Airy function) and cusp caustic (Pearcey function) are
 * classified as Thom catastrophes A2 and A3 respectively. This application
 * computes caustic geometries using the catastrophe theory framework. */
#include "../include/tec_core.h"
#include "../include/tec_bifurcation.h"
#include <math.h>
#include <stdio.h>

/* Airy intensity near fold caustic: I(x) = |Ai(-x)|^2 */
static double airy_intensity(double x) {
    /* Approximate Airy function via catastrophe potential */
    double V = tec_potential_fold(x, -1.0);
    return exp(-0.5 * V * V);  /* simplified approximation */
}

/* Pearcey pattern near cusp caustic */
static double pearcey_intensity(double x, double y) {
    double V = tec_potential_cusp(x, -2.0, y);
    return exp(-0.3 * fabs(V));
}

int tec_optics_caustic(double a_range, double b_range, int n_pts,
                       double* fold_intensity, double* cusp_intensity) {
    int n = 0;

    /* Fold caustic: vary a, compute Airy-like intensity */
    printf("[tec_app1] Fold caustic intensity profile:\n");
    double a_step = 2.0 * a_range / (double)(n_pts - 1);
    for (int i = 0; i < n_pts; i++) {
        double a = -a_range + i * a_step;
        fold_intensity[i] = airy_intensity(a);
        if (i % (n_pts/5) == 0)
            printf("  a=%.3f I=%.6f\n", a, fold_intensity[i]);
    }

    /* Cusp caustic region */
    printf("[tec_app1] Cusp caustic (bifurcation set check):\n");
    TEC_BifurcationSet* cusp_bs = tec_bif_set_cusp(20, a_range);
    for (int i = 0; i < cusp_bs->n_points; i++) {
        double a = cusp_bs->points[i].control[0];
        double b = cusp_bs->points[i].control[1];
        double I = pearcey_intensity(a, b);
        cusp_intensity[i] = I;
        if (i % 4 == 0)
            printf("  a=%.3f b=%.4f I=%.6f\n", a, b, I);
    }
    n = cusp_bs->n_points;
    tec_bif_set_free(cusp_bs);
    return n;
}

int main(void) {
    double f_i[50], c_i[100];
    int n = tec_optics_caustic(3.0, 1.5, 30, f_i, c_i);
    printf("=== tec_app1: Caustics in Optics ===\n");
    printf("Fold profile: max_I=%.4f, Cusp region: %d points\n", f_i[15], n);
    return 0;
}
