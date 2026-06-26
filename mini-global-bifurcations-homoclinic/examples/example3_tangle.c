/* Homoclinic tangle: compute manifold intersections and horseshoe detection */
#include "../include/homoclinic_tangle.h"
#include "../include/homoclinic_orbit.h"
#include <stdio.h>
#include <math.h>

static void pendulum_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double eps = s->params[0], omega = s->params[1];
    dx->x[0] = x.x[1];
    dx->x[1] = -sin(x.x[0]) + eps * sin(omega * s->t);
}
static void pendulum_jac(HOMSystem* s, HOMState x, double* J) {
    J[0] = 0; J[1] = 1;
    J[2] = -cos(x.x[0]); J[3] = 0;
    (void)s;
}

int main(void) {
    printf("=== Homoclinic Tangle: Forced Pendulum ===\n\n");

    HOMSystem* sys = hom_create(2, 2, 0.02);
    sys->rhs = pendulum_rhs;
    sys->jacobian = pendulum_jac;
    hom_set_param(sys, 0, 0.15);   /* eps (forcing amplitude) */
    hom_set_param(sys, 1, 1.0);    /* omega (forcing frequency) */

    double x_saddle[] = {3.141592653589793, 0.0};  /* saddle at (pi, 0) */
    hom_set_state(sys, x_saddle);

    /* Find saddle fixed point */
    HOMFixedPoint fp;
    int fp_rc = hom_find_fixed_point(sys, x_saddle, 1e-8, 50, &fp);
    printf("Saddle FP: found=%s, type=%d\n", fp_rc >= 0 ? "YES" : "NO", fp.type);

    /* Compute stable and unstable manifolds */
    HOMTrajectory Ws, Wu;
    int n_ws = hom_compute_stable_manifold(sys, &fp, 0.01, 8.0, 500, &Ws);
    int n_wu = hom_compute_unstable_manifold(sys, &fp, 0.01, 8.0, 500, &Wu);
    printf("Manifolds: Ws=%d pts, Wu=%d pts\n", n_ws, n_wu);

    /* Detect manifold intersection */
    HOMState inter;
    bool intersect = hom_manifolds_intersect(&Ws, &Wu, 0.05, &inter);
    printf("Manifolds intersect: %s", intersect ? "YES" : "NO");
    if (intersect) printf(" at (%.4f, %.4f)", inter.x[0], inter.x[1]);
    printf("\n");

    /* Tangle analysis */
    HOMTangle* tangle = hom_tangle_create();
    int td = hom_tangle_detect(sys, &fp, tangle, 0.01, 15.0, 300);
    printf("Tangle: %d intersections, has_tangle=%s\n",
           tangle->n_intersections, tangle->has_tangle ? "YES" : "NO");

    /* Horseshoe detection */
    if (tangle->has_tangle) {
        bool hs = hom_tangle_has_horseshoe(tangle);
        printf("Smale horseshoe: %s\n", hs ? "DETECTED" : "not found");
        if (hs) {
            printf("Horseshoe count: %d\n", tangle->smale_horseshoe_count);
            bool chaos = hom_horseshoe_implies_chaos(tangle);
            printf("Horseshoe => chaos: %s\n", chaos ? "YES (Smale-Birkhoff)" : "NO");
        }
        double fd = hom_tangle_fractal_dimension(tangle);
        printf("Fractal dimension: %.4f\n", fd);
    }

    /* Lobe analysis */
    HOMLobe lobe;
    bool has_lobe = hom_tangle_get_lobe(tangle, 0, &lobe);
    if (has_lobe) {
        printf("Lobe 0: center=(%.3f,%.3f) area=%.4f chaotic=%s\n",
               lobe.center.x[0], lobe.center.x[1], lobe.area,
               lobe.is_chaotic ? "YES" : "NO");
    }

    hom_tangle_print(tangle);
    hom_tangle_free(tangle);
    hom_free(sys);
    printf("\nExample 3 PASSED\n");
    return 0;
}
