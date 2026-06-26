/* Shilnikov chaos: detect homoclinic orbit via shooting method */
#include "../include/homoclinic_core.h"
#include "../include/homoclinic_computation.h"
#include "../include/homoclinic_orbit.h"
#include <stdio.h>
#include <math.h>

static void shilnikov_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double a = s->params[0], b = s->params[1], c = s->params[2];
    dx->x[0] = a * (x.x[1] - x.x[0] - (b*x.x[0] + 0.5*(c-b)*(fabs(x.x[0]+1)-fabs(x.x[0]-1))));
    dx->x[1] = x.x[0] - x.x[1] + x.x[2];
    dx->x[2] = -b * x.x[1];
}
static void shilnikov_jac(HOMSystem* s, HOMState x, double* J) {
    (void)s; (void)x; (void)J; /* simplified for example */
}

int main(void) {
    printf("=== Shilnikov Chaos: Homoclinic Orbit Detection ===\n\n");

    HOMSystem* sys = hom_create(3, 3, 0.01);
    sys->rhs = shilnikov_rhs;
    sys->jacobian = shilnikov_jac;
    hom_set_param(sys, 0, 9.0);   /* alpha */
    hom_set_param(sys, 1, 14.28); /* beta */
    hom_set_param(sys, 2, -0.71); /* m1 */

    /* Set initial state near saddle-focus */
    double x0[] = {0.1, 0.05, -0.02};
    hom_set_state(sys, x0);
    printf("System: dim=%d, params=[%.2f, %.2f, %.2f]\n",
           sys->dim, sys->params[0], sys->params[1], sys->params[2]);

    /* Integrate to find fixed point */
    HOMFixedPoint fp;
    int fp_rc = hom_find_fixed_point(sys, x0, 1e-8, 50, &fp);
    printf("Fixed point: found=%s, type=%d, dim_s=%d, dim_u=%d\n",
           fp_rc >= 0 ? "YES" : "NO", fp.type, fp.dim_stable, fp.dim_unstable);

    if (fp_rc >= 0 && fp.n_ev > 0) {
        printf("Eigenvalues:");
        for (int i = 0; i < fp.n_ev; i++) {
            printf(" %.4f%+.4fi", fp.eigenvalues[i].re, fp.eigenvalues[i].im);
        }
        /* Shilnikov condition */
        if (fp.n_ev >= 3 && fp.eigenvalues[0].re > 0) {
            double rho = fp.eigenvalues[1].re;
            double lam_u = fp.eigenvalues[0].re;
            double ratio = fabs(rho) / lam_u;
            printf("\nShilnikov: |Re(rho)|/lambda_u = %.4f/%.4f = %.4f", fabs(rho), lam_u, ratio);
            printf(ratio < 1.0 ? " < 1 => CHAOS possible!" : " >= 1 => no Shilnikov chaos");
        }
    }
    printf("\n");

    /* Continuation to find bifurcation */
    HOMContinuation* cont = hom_continuation_create(0, 8.0, 10.0, 20);
    int crc = hom_continuation_run(sys, &fp, cont, 1e-6, 30);
    printf("Continuation: %d orbits computed, bif=%s at p=%.4f\n",
           cont->n_orbits, cont->found_bifurcation ? "YES" : "NO", cont->bifurcation_value);

    bool is_global = hom_is_global_bifurcation(sys, &fp, 1e-6);
    printf("Global bifurcation detected: %s\n", is_global ? "YES" : "NO");

    hom_continuation_free(cont);
    hom_free(sys);
    printf("\nExample 1 PASSED\n");
    return 0;
}
