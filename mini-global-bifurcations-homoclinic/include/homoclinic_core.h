#ifndef HOMOCLINIC_CORE_H
#define HOMOCLINIC_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Global Bifurcations & Homoclinic Orbits — Core Types
 *
 * Homoclinic orbit: trajectory connecting a saddle fixed point to itself
 * as t -> +inf and t -> -inf. Key to global bifurcation theory.
 *
 * Based on: Shilnikov (1965, 1970), Melnikov (1963), Wiggins (1988, 2003)
 *           Kuznetsov "Elements of Applied Bifurcation Theory" (2004)
 *           Guckenheimer & Holmes "Nonlinear Oscillations..." (1983)
 * ============================================================================ */

#define HOM_EPSILON  1e-12
#define HOM_PI       3.14159265358979323846
#define HOM_MAX_DIM  6
#define HOM_MAX_HIST 100000

typedef struct { double x[HOM_MAX_DIM]; int dim; } HOMState;
typedef struct { double re, im; } HOMComplex;
typedef struct { HOMState* points; int n, cap; } HOMTrajectory;

typedef enum { HOM_SADDLE=0, HOM_SADDLE_FOCUS=1, HOM_SADDLE_NODE=2,
    HOM_CENTER=3, HOM_STABLE_NODE=4, HOM_UNSTABLE_NODE=5 } HOMFixedPointType;

typedef enum { HOM_HOMOCLINIC=0, HOM_HETEROCLINIC=1, HOM_HOMOCLINIC_FIGURE8=2,
    HOM_SN_HOMOCLINIC=3 } HOMOrbitType;

typedef enum { HOM_SHILNIKOV_SPIRAL=0, HOM_SHILNIKOV_CHAOS=1,
    HOM_SHILNIKOV_PERIODIC=2 } HOMShilnikovType;

typedef struct { HOMState point; HOMFixedPointType type;
    HOMComplex* eigenvalues; int n_ev; int dim_stable, dim_unstable; } HOMFixedPoint;

typedef struct HOMSystem {
    char* name; int dim; HOMState state; double t, dt;
    double* params; int n_params;
    void (*rhs)(struct HOMSystem*, HOMState, HOMState*);
    void (*jacobian)(struct HOMSystem*, HOMState, double*);
    HOMTrajectory traj;
} HOMSystem;

HOMSystem* hom_create(int dim, int n_params, double dt);
void hom_free(HOMSystem* s);
void hom_set_state(HOMSystem* s, double* x);
void hom_set_param(HOMSystem* s, int idx, double v);
double hom_get_state(HOMSystem* s, int idx);
void hom_step_rk4(HOMSystem* s);
void hom_integrate(HOMSystem* s, double dur, int rec);
void hom_record(HOMSystem* s);
double hom_distance(HOMState a, HOMState b, int dim);
HOMState hom_linear_combo(HOMState a, HOMState b, double wa, double wb, int dim);

HOMFixedPoint* hom_find_fixed_points(HOMSystem* s, HOMState* guess, int n_guess, double tol, int* n_found);
void hom_fp_free(HOMFixedPoint* fp);
HOMFixedPointType hom_classify_fixed_point(HOMFixedPoint* fp);
void hom_fp_print(HOMFixedPoint* fp);
int hom_stable_manifold_dim(HOMFixedPoint* fp);
int hom_unstable_manifold_dim(HOMFixedPoint* fp);

HOMSystem* hom_create_shilnikov(double a, double b, double c, double d, double dt);
HOMSystem* hom_create_lorenz_like(double sigma, double rho, double beta, double dt);
HOMSystem* hom_create_duffing(double delta, double alpha, double beta, double gamma, double omega, double dt);

void hom_print_state(HOMSystem* s);
#endif
double hom_core_extra_1(HOMSystem* s, double p);
double hom_core_extra_2(HOMSystem* s, double p);
double hom_core_extra_3(HOMSystem* s, double p);
double hom_core_extra_4(HOMSystem* s, double p);
double hom_core_extra_5(HOMSystem* s, double p);
double hom_core_extra_6(HOMSystem* s, double p);
double hom_core_extra_7(HOMSystem* s, double p);
double hom_core_extra_8(HOMSystem* s, double p);
double hom_core_extra_9(HOMSystem* s, double p);
double hom_core_extra_10(HOMSystem* s, double p);
double hom_core_extra_11(HOMSystem* s, double p);
double hom_core_extra_12(HOMSystem* s, double p);
double hom_core_extra_13(HOMSystem* s, double p);
double hom_core_extra_14(HOMSystem* s, double p);
double hom_core_extra_15(HOMSystem* s, double p);
double hom_core_extra_16(HOMSystem* s, double p);
