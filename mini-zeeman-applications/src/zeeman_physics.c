#include "zeeman_physics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Euler Buckling (Zeeman 1976)
 *
 * An elastic beam under axial compression P buckles when P
 * exceeds the critical (Euler) load P_cr.
 *
 * For a pinned-pinned beam: P_cr = pi^2 * EI / L^2
 * For a clamped-clamped beam: P_cr = 4 * pi^2 * EI / L^2
 *
 * The buckling is a pitchfork bifurcation (Z_2 symmetric cusp):
 *
 * Potential energy: V(w) = (EI/2)*int(w'')^2 dx - (P/2)*int(w')^2 dx
 *
 * Leading to amplitude equation:
 *   dA/dt = (P/P_cr - 1)*A - (something)*A^3
 *
 * This is the normal form of a supercritical pitchfork bifurcation,
 * equivalent to the cusp catastrophe with b=0.
 *
 * Zeeman showed that Euler buckling is the simplest physical
 * demonstration of a catastrophe.
 * ============================================================== */

EulerBuckling* euler_buckling_create(double load, double length, double EI,
                                      int boundary) {
    EulerBuckling* e = calloc(1, sizeof(EulerBuckling));
    if (!e) return NULL;
    e->load = load;
    e->deflection = 0.01;   /* small initial imperfection */
    e->length = length;
    e->EI = EI;
    e->boundary_cond = (double)boundary;
    e->t = 0.0;
    return e;
}

void euler_buckling_free(EulerBuckling* e) {
    free(e);
}

double euler_buckling_critical_load(const EulerBuckling* e) {
    /* Euler buckling formula */
    double factor = (e->boundary_cond < 0.5) ? 1.0 : 4.0;
    return factor * ZM_PI * ZM_PI * e->EI / (e->length * e->length);
}

double euler_buckling_potential(const EulerBuckling* e, double deflection) {
    /* Post-buckling potential: V(A) ~ (P_cr - P)*A^2/2 + (EI)*A^4/4
     * This is exactly the cusp form with:
     * a = (P_cr - P) / P_cr  (normal factor)
     * b = 0 (symmetric) */
    double Pcr = euler_buckling_critical_load(e);
    double a = (Pcr - e->load) / Pcr;
    double b = 0.0;
    return cusp_potential(deflection, a, b) * e->EI / e->length;
}

double euler_buckling_deflection(const EulerBuckling* e, double load) {
    /* Post-buckling deflection amplitude (Koiter, 1945)
     * For P > Pcr: A = sqrt(8*(P-Pcr)/(3*Pcr)) * L/pi */
    double Pcr = euler_buckling_critical_load(e);
    if (load <= Pcr) return 0.0;
    double ratio = (load - Pcr) / Pcr;
    return sqrt(8.0 * ratio / 3.0) * e->length / ZM_PI;
}

void euler_buckling_step(EulerBuckling* e, double dt) {
    double Pcr = euler_buckling_critical_load(e);
    double a = (Pcr - e->load) / Pcr;

    /* Dynamics: dA/dt = -dV/dA = -(A^3 + a*A) */
    double force = -(e->deflection * e->deflection * e->deflection
                    + a * e->deflection);
    e->deflection += force * dt;
    e->t += dt;

    if (e->deflection > 2.0 * e->length) e->deflection = 2.0 * e->length;
    if (e->deflection < -2.0 * e->length) e->deflection = -2.0 * e->length;
}

void euler_buckling_simulate(EulerBuckling* e, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        euler_buckling_step(e, dt);
    }
}

double euler_buckling_post_buckling_stiffness(const EulerBuckling* e) {
    /* Post-buckling tangent stiffness */
    double Pcr = euler_buckling_critical_load(e);
    if (e->load <= Pcr) return e->EI;
    double ratio = (e->load - Pcr) / Pcr;
    return e->EI / (1.0 + 3.0 * ratio);
}

bool euler_buckling_has_buckled(const EulerBuckling* e) {
    double Pcr = euler_buckling_critical_load(e);
    return e->load > Pcr && fabs(e->deflection) > e->length * 1e-4;
}

void euler_buckling_print(const EulerBuckling* e) {
    double Pcr = euler_buckling_critical_load(e);
    printf("EulerBuckling: P=%.2f Pcr=%.2f defl=%.6f L=%.2f EI=%.2f t=%.2f\n",
           e->load, Pcr, e->deflection, e->length, e->EI, e->t);
    printf("  Load ratio: %.2f  Buckled: %s  Stiffness: %.4f\n",
           e->load / Pcr,
           euler_buckling_has_buckled(e) ? "YES" : "NO",
           euler_buckling_post_buckling_stiffness(e));
}

/* ==============================================================
 * Phase Transition (van der Waals gas)
 *
 * The van der Waals equation of state:
 *   (P + a/V^2)(V - b) = RT
 *
 * rearranged: P(V,T) = RT/(V-b) - a/V^2
 *
 * Near the critical point (Tc, Pc, Vc), Taylor expansion yields
 * a cusp catastrophe form.
 *
 * Critical parameters:
 *   Vc = 3b
 *   Tc = 8a/(27Rb)
 *   Pc = a/(27b^2)
 *
 * The Helmholtz free energy F(V,T) serves as the potential.
 * Below Tc, F has two minima (liquid and vapor) and one maximum.
 *
 * Maxwell construction (equal-area rule) determines the
 * coexistence pressure where liquid and vapor have equal
 * chemical potential.
 * ============================================================== */

PhaseTransition* phase_transition_create(double T, double V, double a_vdw,
                                          double b_vdw) {
    PhaseTransition* p = calloc(1, sizeof(PhaseTransition));
    if (!p) return NULL;
    p->T = T;
    p->V = V;
    p->a_vdw = a_vdw;
    p->b_vdw = b_vdw;
    p->R_gas = 8.314462618;
    p->t = 0.0;
    phase_transition_compute_critical(p);
    /* Initial pressure from EOS */
    p->P = phase_transition_pressure(p, V);
    return p;
}

void phase_transition_free(PhaseTransition* p) {
    free(p);
}

double phase_transition_pressure(const PhaseTransition* p, double V) {
    if (V <= p->b_vdw) return 1e10;  /* avoid singularity */
    return p->R_gas * p->T / (V - p->b_vdw) - p->a_vdw / (V * V);
}

double phase_transition_potential(const PhaseTransition* p, double V) {
    /* Helmholtz free energy relative to ideal gas:
     * F(V,T) = -RT*ln(V-b) - a/V + F_id(T) */
    if (V <= p->b_vdw) return 1e10;
    return -p->R_gas * p->T * log(V - p->b_vdw) - p->a_vdw / V;
}

void phase_transition_isotherm(PhaseTransition* p, double* V_vals,
                                double* P_vals, int n) {
    double V_min = p->b_vdw * 1.1;
    double V_max = 10.0 * p->Vc;
    for (int i = 0; i < n; i++) {
        double V = V_min + (V_max - V_min) * i / (n - 1);
        V_vals[i] = V;
        P_vals[i] = phase_transition_pressure(p, V);
    }
}

void phase_transition_compute_critical(PhaseTransition* p) {
    /* Critical point from van der Waals parameters */
    p->Vc = 3.0 * p->b_vdw;
    p->Tc = 8.0 * p->a_vdw / (27.0 * p->R_gas * p->b_vdw);
    p->Pc = p->a_vdw / (27.0 * p->b_vdw * p->b_vdw);
}

double phase_transition_maxwell_construction(const PhaseTransition* p,
                                              double* Vl, double* Vg) {
    /* Maxwell equal-area rule:
     * Integral[P(V) - P_eq, {V, Vl, Vg}] = 0
     * P_eq * (Vg - Vl) = Integral[P(V), {V, Vl, Vg}]
     *
     * For van der Waals: analytical solution via
     * P_eq computed iteratively. */
    if (p->T >= p->Tc) {
        if (Vl) *Vl = p->Vc;
        if (Vg) *Vg = p->Vc;
        return p->Pc;
    }

    /* Bisection search for coexistence pressure */
    double P_low = 0.0;
    double P_high = phase_transition_pressure(p, p->Vc);
    double P_eq = (P_low + P_high) / 2.0;

    for (int iter = 0; iter < 50; iter++) {
        /* Find volumes at this pressure */
        /* Solve P(V) = P_eq -> cubic in V */
        /* (P + a/V^2)(V - b) = RT
         * P*V^3 - (P*b + RT)*V^2 + a*V - a*b = 0 */
        double c3 = P_eq;
        double c2 = -(P_eq * p->b_vdw + p->R_gas * p->T);
        double c1 = p->a_vdw;
        double c0 = -p->a_vdw * p->b_vdw;

        /* Solve cubic: c3*V^3 + c2*V^2 + c1*V + c0 = 0 */
        /* Normalize to depressed cubic */
        double b2 = c2 / c3;
        double b1 = c1 / c3;
        double b0 = c0 / c3;

        /* Depress: let V = y - b2/3 */
        double p_dep = b1 - b2 * b2 / 3.0;
        double q_dep = b0 - b1 * b2 / 3.0 + 2.0 * b2 * b2 * b2 / 27.0;

        double roots[3];
        int n = zeeman_solve_cubic_cardano(p_dep, q_dep, roots);
        /* Shift back */
        double shift = b2 / 3.0;
        for (int i = 0; i < n; i++) roots[i] -= shift;

        /* Sort and find reasonable roots */
        if (n >= 2) {
            double V1 = roots[0], V3 = roots[n-1];
            /* Check equal-area condition */
            double area = phase_transition_pressure(p, (V1+V3)/2)
                        * (V3 - V1);
            double integral = 0.0;
            int ni = 100;
            double dV = (V3 - V1) / ni;
            for (int i = 1; i < ni; i++) {
                double V = V1 + i * dV;
                integral += phase_transition_pressure(p, V) * dV;
            }
            if (integral > area) P_high = P_eq;
            else P_low = P_eq;
        }
        P_eq = (P_low + P_high) / 2.0;
    }

    /* Final root finding at P_eq */
    double c3 = P_eq;
    double c2 = -(P_eq * p->b_vdw + p->R_gas * p->T);
    double c1 = p->a_vdw;
    double c0 = -p->a_vdw * p->b_vdw;
    double b2 = c2 / c3, b1 = c1 / c3, b0 = c0 / c3;
    double pd = b1 - b2*b2/3.0, qd = b0 - b1*b2/3.0 + 2.0*b2*b2*b2/27.0;
    double roots[3];
    int n = zeeman_solve_cubic_cardano(pd, qd, roots);
    double shift = b2 / 3.0;
    for (int i = 0; i < n; i++) roots[i] -= shift;

    if (n >= 3) {
        if (Vl) *Vl = roots[0];  /* liquid volume (smallest) */
        if (Vg) *Vg = roots[2];  /* vapor volume (largest) */
    } else if (Vl && Vg) {
        *Vl = *Vg = p->Vc;
    }
    return P_eq;
}

void phase_transition_spinodal(const PhaseTransition* p,
                                double* V_sp1, double* V_sp2) {
    /* Spinodal: dP/dV = 0
     * -RT/(V-b)^2 + 2a/V^3 = 0
     * 2a*(V-b)^2 = RT*V^3
     * This is a cubic in V. */
    if (p->T >= p->Tc) {
        if (V_sp1) *V_sp1 = p->Vc;
        if (V_sp2) *V_sp2 = p->Vc;
        return;
    }

    /* Solve numerically by scanning */
    double V_min = p->b_vdw * 1.01;
    double V_max = p->Vc * 5.0;
    int found = 0;
    double prev_dP = 0.0;
    for (int i = 1; i < 1000; i++) {
        double V = V_min + (V_max - V_min) * i / 999.0;
        double V_eps = V * 1.001;
        double dP = (phase_transition_pressure(p, V_eps)
                   - phase_transition_pressure(p, V)) / (V_eps - V);
        if (i > 1 && prev_dP * dP <= 0.0 && fabs(dP) < 100.0) {
            if (found == 0 && V_sp1) *V_sp1 = V;
            else if (found == 1 && V_sp2) *V_sp2 = V;
            found++;
        }
        prev_dP = dP;
    }
    if (found < 2 && V_sp2) *V_sp2 = *V_sp1;
}

void phase_transition_step(PhaseTransition* p, double dt) {
    /* V relaxes to equilibrium at constant T */
    if (p->T >= p->Tc) {
        /* Supercritical: single phase */
        p->V += (p->Vc - p->V) * 0.1 * dt;
    } else {
        /* Subcritical: move toward nearest stable phase */
        double Vl, Vg;
        phase_transition_maxwell_construction(p, &Vl, &Vg);
        if (p->V < (Vl + Vg) / 2.0) {
            p->V += (Vl - p->V) * 0.1 * dt;
        } else {
            p->V += (Vg - p->V) * 0.1 * dt;
        }
    }
    p->P = phase_transition_pressure(p, p->V);
    p->t += dt;
}

void phase_transition_simulate(PhaseTransition* p, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        phase_transition_step(p, dt);
    }
}

bool phase_transition_is_supercritical(const PhaseTransition* p) {
    return p->T >= p->Tc;
}

void phase_transition_print(const PhaseTransition* p) {
    printf("PhaseTransition: T=%.2f P=%.4f V=%.4f a=%.2f b=%.2f t=%.2f\n",
           p->T, p->P, p->V, p->a_vdw, p->b_vdw, p->t);
    printf("  Critical: Tc=%.2f Pc=%.4f Vc=%.4f\n", p->Tc, p->Pc, p->Vc);
    printf("  Regime: %s (T/Tc=%.3f)\n",
           phase_transition_is_supercritical(p) ? "SUPERCRITICAL" : "SUB-CRITICAL",
           p->T / p->Tc);
    if (!phase_transition_is_supercritical(p)) {
        double Vl, Vg, Peq;
        Peq = phase_transition_maxwell_construction(p, &Vl, &Vg);
        printf("  Coexistence: Pl=%.4f Pg=%.4f Vl=%.4f Vg=%.4f Peq=%.4f\n",
               phase_transition_pressure(p, Vl),
               phase_transition_pressure(p, Vg), Vl, Vg, Peq);
    }
}

/* ==============================================================
 * Ship Stability (Zeeman 1977)
 *
 * Ship roll dynamics exhibit a cusp catastrophe at the point
 * of vanishing metacentric height (GM).
 *
 * Righting moment: M(phi) = Delta * GM * sin(phi)
 *   where Delta = displacement, GM = metacentric height
 *
 * When GM > 0: stable equilibrium at phi = 0
 * When GM < 0: ship capsizes or finds new equilibrium at large angle
 * When GM = 0: cusp catastrophe point
 *
 * The potential energy: V(phi) = Delta * GM * (1 - cos(phi))
 * Plus nonlinear terms at large angles.
 *
 * Nonlinear righting arm: GZ(phi) = GM * sin(phi)
 *   + (BM/2)*tan^2(phi)*sin(phi) (wall-sided formula)
 *
 * The cusp catastrophe appears when considering the effect
 * of beam winds or waves as a control parameter b.
 * ============================================================== */

ShipStability* ship_stability_create(double gm, double beam, double damping) {
    ShipStability* s = calloc(1, sizeof(ShipStability));
    if (!s) return NULL;
    s->roll_angle = 0.01;       /* small initial heel */
    s->roll_rate = 0.0;
    s->metacentric_height = gm;
    s->beam = beam;
    s->damping_coef = damping;
    s->wave_amplitude = 0.0;
    s->wave_frequency = 0.5;
    s->t = 0.0;
    return s;
}

void ship_stability_free(ShipStability* s) {
    free(s);
}

double ship_stability_righting_moment(const ShipStability* s, double angle) {
    /* GZ(phi) = GM*sin(phi) + nonlinear correction
     * Wall-sided formula: GZ(phi) = sin(phi)*(GM + BM*tan^2(phi)/2)
     * where BM is approximated from beam */
    double BM = s->beam * s->beam / 12.0;  /* approximate BM for box-shaped hull */
    double tan_phi = tan(angle);
    double sin_phi = sin(angle);

    /* Limit tan for large angles */
    if (fabs(tan_phi) > 10.0) tan_phi = (tan_phi > 0) ? 10.0 : -10.0;

    double GZ = sin_phi * (s->metacentric_height + 0.5 * BM * tan_phi * tan_phi);
    return GZ;  /* normalized by displacement */
}

double ship_stability_potential(const ShipStability* s, double angle) {
    /* V(phi) = -Integral[GZ(phi') dphi', 0, phi]
     * For small angles: V(phi) = GM*(1 - cos(phi))
     * Full potential from integration */
    double V = 0.0;
    int n_int = 100;
    double dphi = angle / n_int;
    for (int i = 0; i < n_int; i++) {
        double phi_mid = (i + 0.5) * dphi;
        V -= ship_stability_righting_moment(s, phi_mid) * dphi;
    }
    return V;
}

CuspEquilibrium* ship_stability_equilibrium(const ShipStability* s) {
    /* Find equilibrium roll angles where GZ(phi) = 0 */
    /* For GM > 0: phi = 0 is the primary stable equilibrium */
    /* We'll treat it as cusp: a = -GM, b = wind/wave moment */
    double a = -s->metacentric_height * 10.0;  /* scale to cusp range */
    double b = s->wave_amplitude * 5.0;
    return cusp_equilibrium(a, b);
}

void ship_stability_step(ShipStability* s, double dt) {
    /* Roll equation: I*phi'' + B*phi' + Delta*GZ(phi) = M_wave(t) */
    double inertia = 1.0;
    double GZ = ship_stability_righting_moment(s, s->roll_angle);
    double wave_moment = s->wave_amplitude * sin(s->wave_frequency * s->t);
    double damping_torque = s->damping_coef * s->roll_rate;

    double alpha = (wave_moment - damping_torque - GZ) / inertia;
    s->roll_rate += alpha * dt;
    s->roll_angle += s->roll_rate * dt;
    s->t += dt;

    /* Keep angle reasonable */
    while (s->roll_angle > ZM_PI) s->roll_angle -= ZM_TWO_PI;
    while (s->roll_angle < -ZM_PI) s->roll_angle += ZM_TWO_PI;
}

void ship_stability_simulate(ShipStability* s, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        ship_stability_step(s, dt);
    }
}

bool ship_stability_capsize_check(const ShipStability* s, double max_angle) {
    /* Capsize if roll angle exceeds max_angle or GM negative */
    return fabs(s->roll_angle) > max_angle || s->metacentric_height <= 0.0;
}

double ship_stability_critical_gm(const ShipStability* s) {
    /* Critical GM where stability vanishes:
     * For a given wave/wind moment, solve for minimum GM.
     * Approximate: GM_crit ~ wave_amplitude / sin(max_angle) */
    double max_safe_angle = 0.5;  /* ~30 degrees */
    return s->wave_amplitude / sin(max_safe_angle);
}

void ship_stability_set_wave(ShipStability* s, double amplitude,
                              double frequency) {
    s->wave_amplitude = amplitude;
    s->wave_frequency = frequency;
}

void ship_stability_print(const ShipStability* s) {
    printf("ShipStability: roll=%.2f deg rate=%.3f GM=%.4f beam=%.2f damp=%.2f t=%.2f\n",
           s->roll_angle * ZM_RAD2DEG, s->roll_rate,
           s->metacentric_height, s->beam, s->damping_coef, s->t);
    double GZ = ship_stability_righting_moment(s, s->roll_angle);
    printf("  GZ=%.6f  CapsizeRisk: %s  CriticalGM: %.4f\n",
           GZ,
           ship_stability_capsize_check(s, 1.0) ? "HIGH" : "LOW",
           ship_stability_critical_gm(s));
}