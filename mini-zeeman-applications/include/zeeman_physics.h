#ifndef ZEEMAN_PHYSICS_H
#define ZEEMAN_PHYSICS_H
#include "zeeman_core.h"

/* Zeeman's Physics & Engineering Applications */

typedef struct {
    double load;            /* axial compressive load */
    double deflection;      /* lateral deflection amplitude */
    double length;          /* beam length */
    double EI;              /* flexural rigidity */
    double boundary_cond;   /* 0=pinned-pinned, 1=clamped-clamped */
    double t;
} EulerBuckling;

typedef struct {
    double T;           /* temperature */
    double P;           /* pressure */
    double V;           /* specific volume */
    double a_vdw;       /* van der Waals attraction */
    double b_vdw;       /* excluded volume */
    double R_gas;       /* gas constant */
    double Tc, Pc, Vc;  /* critical point */
    double t;
} PhaseTransition;

typedef struct {
    double roll_angle;          /* roll angle (radians) */
    double roll_rate;           /* angular velocity */
    double metacentric_height;  /* GM: stability parameter */
    double beam;                /* ship beam (width) */
    double damping_coef;        /* roll damping */
    double wave_amplitude;
    double wave_frequency;
    double t;
} ShipStability;

/* Euler buckling API (Zeeman 1976) */
EulerBuckling* euler_buckling_create(double load, double length, double EI, int boundary);
void           euler_buckling_free(EulerBuckling* e);
double         euler_buckling_critical_load(const EulerBuckling* e);
double         euler_buckling_potential(const EulerBuckling* e, double deflection);
double         euler_buckling_deflection(const EulerBuckling* e, double load);
void           euler_buckling_step(EulerBuckling* e, double dt);
void           euler_buckling_simulate(EulerBuckling* e, double duration, double dt);
double         euler_buckling_post_buckling_stiffness(const EulerBuckling* e);
bool           euler_buckling_has_buckled(const EulerBuckling* e);
void           euler_buckling_print(const EulerBuckling* e);

/* Phase transition API */
PhaseTransition* phase_transition_create(double T, double V, double a_vdw, double b_vdw);
void             phase_transition_free(PhaseTransition* p);
double           phase_transition_pressure(const PhaseTransition* p, double V);
double           phase_transition_potential(const PhaseTransition* p, double V);
void             phase_transition_isotherm(PhaseTransition* p, double* V_vals, double* P_vals, int n);
void             phase_transition_compute_critical(PhaseTransition* p);
double           phase_transition_maxwell_construction(const PhaseTransition* p, double* Vl, double* Vg);
void             phase_transition_spinodal(const PhaseTransition* p, double* V_sp1, double* V_sp2);
void             phase_transition_step(PhaseTransition* p, double dt);
void             phase_transition_simulate(PhaseTransition* p, double duration, double dt);
bool             phase_transition_is_supercritical(const PhaseTransition* p);
void             phase_transition_print(const PhaseTransition* p);

/* Ship stability API (Zeeman 1977) */
ShipStability* ship_stability_create(double gm, double beam, double damping);
void           ship_stability_free(ShipStability* s);
double         ship_stability_righting_moment(const ShipStability* s, double angle);
double         ship_stability_potential(const ShipStability* s, double angle);
CuspEquilibrium* ship_stability_equilibrium(const ShipStability* s);
void           ship_stability_step(ShipStability* s, double dt);
void           ship_stability_simulate(ShipStability* s, double duration, double dt);
bool           ship_stability_capsize_check(const ShipStability* s, double max_angle);
double         ship_stability_critical_gm(const ShipStability* s);
void           ship_stability_set_wave(ShipStability* s, double amplitude, double frequency);
void           ship_stability_print(const ShipStability* s);

#endif /* ZEEMAN_PHYSICS_H */