/* hom_app2.c -- L7-2: Power grid transient stability (Chiang 1988, BCU method).
 * A single-machine-infinite-bus (SMIB) power system:
 * M*delta_ddot + D*delta_dot + P_max*(sin(delta) - sin(delta_s)) = 0
 * has a homoclinic orbit at the critical clearing angle delta_cr.
 * Beyond this boundary, the system escapes the potential well
 * (loss of synchronism = blackout). */
#include "../include/homoclinic_core.h"
#include "../include/homoclinic_orbit.h"
#include <math.h>
#include <stdio.h>

static void smib_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double M = s->params[0], D = s->params[1], P_max = s->params[2];
    double delta_s = s->params[3];
    dx->x[0] = x.x[1];  /* d(delta)/dt = omega */
    dx->x[1] = (P_max * (sin(delta_s) - sin(x.x[0])) - D * x.x[1]) / M;
}
static void smib_jac(HOMSystem* s, HOMState x, double* J) {
    double M = s->params[0], D = s->params[1], P_max = s->params[2];
    J[0] = 0.0; J[1] = 1.0;
    J[2] = -P_max * cos(x.x[0]) / M; J[3] = -D / M;
}

int hom_power_grid(double P_max, double V_inf, double X_line,
                   double* crit_clearing_angle, double* crit_clearing_time,
                   double* stability_margin) {
    double M = 0.05;       /* inertia constant [pu] */
    double D = 0.02;       /* damping */
    double P_mech = 1.0;   /* mechanical power */
    double delta_s = asin(P_mech * X_line / (V_inf * V_inf));
    double P_elec_max = V_inf * V_inf / X_line;

    HOMSystem* sys = hom_create(2, 4, 0.005);
    sys->rhs = smib_rhs;
    sys->jacobian = smib_jac;
    hom_set_param(sys, 0, M);
    hom_set_param(sys, 1, D);
    hom_set_param(sys, 2, P_elec_max);
    hom_set_param(sys, 3, delta_s);

    /* Saddle point at (pi - delta_s, 0) */
    double x_saddle[] = {M_PI - delta_s, 0.0};
    hom_set_state(sys, x_saddle);

    /* Find fixed point and analyze */
    HOMFixedPoint fp;
    int fp_rc = hom_find_fixed_point(sys, x_saddle, 1e-8, 50, &fp);

    /* Critical clearing angle: energy at saddle = energy at clearing */
    double delta_u = M_PI - delta_s;
    *crit_clearing_angle = acos((delta_u - delta_s) * sin(delta_s) + cos(delta_u));

    /* Critical clearing time via equal area criterion */
    double omega0 = sqrt(P_elec_max * cos(delta_s) / M);
    *crit_clearing_time = sqrt(2.0 * M * (*crit_clearing_angle - delta_s) /
                                (P_mech * omega0));

    /* Stability margin: distance from stable equilibrium to saddle */
    *stability_margin = delta_u - delta_s;

    printf("[hom_app2] Power grid SMIB: V=%.2f X=%.2f\n", V_inf, X_line);
    printf("  Stable eq: delta_s=%.4f rad (%.1f deg)\n", delta_s, delta_s * 180.0 / M_PI);
    printf("  Saddle: delta_u=%.4f rad (%.1f deg)\n", delta_u, delta_u * 180.0 / M_PI);
    printf("  Critical clearing: angle=%.4f rad (%.1f deg), time=%.4f s\n",
           *crit_clearing_angle, *crit_clearing_angle * 180.0 / M_PI,
           *crit_clearing_time);
    printf("  Stability margin: %.4f rad\n", *stability_margin);

    /* Check if near global bifurcation */
    HOMOrbit* orb = hom_orbit_create(500);
    double guess[] = {M_PI - delta_s - 0.1, 0.05};
    HOMState gs = {.dim = 2};
    gs.x[0] = guess[0]; gs.x[1] = guess[1];
    int orb_rc = hom_shoot_homoclinic(sys, &fp, &gs, 1, 1e-6, 80, orb);
    bool is_global = hom_is_global_bifurcation(sys, &fp, 1e-6);
    printf("  Homoclinic orbit found: %s, global bif: %s\n",
           orb_rc >= 0 ? "YES" : "NO", is_global ? "YES" : "NO");

    hom_orbit_free(orb);
    hom_free(sys);
    return 0;
}

int main(void) {
    double d_cr, t_cr, margin;
    hom_power_grid(2.5, 1.0, 0.4, &d_cr, &t_cr, &margin);
    printf("=== hom_app2: Power Grid Transient Stability ===\n");
    printf("delta_cr=%.4f t_cr=%.4fs\n", d_cr, t_cr);
    return 0;
}
