#ifndef ZEEMAN_DYNAMICS_H
#define ZEEMAN_DYNAMICS_H
#include "zeeman_core.h"

/* Zeeman Dynamics - Bifurcation Analysis, Hysteresis, Jumps */

typedef struct {
    double *x_final;
    int     n_trajectories;
    double  divergence_measure;
    double  sensitivity;
} DivergenceTest;

typedef struct {
    double x;       /* fast variable */
    double y;       /* slow variable */
    double epsilon; /* time-scale ratio */
    double a0, b0;  /* baseline cusp parameters */
    double alpha;   /* coupling: a = a0 + alpha*y */
    double beta;    /* slow dynamics parameter */
    double gamma;   /* slow dynamics parameter */
    double t;
} SlowFastCusp;

/* Jump detection */
double zeeman_find_delay_jump(double a_start, double a_end, double b_fixed, int steps, double x0);
double zeeman_find_maxwell_jump(double a_start, double a_end, double b_fixed, int steps);
bool   zeeman_detect_bifurcation_crossing(const double* a, const double* b, int n, int* crossing_idx);
void   zeeman_equilibrium_manifold_section(double a_min, double a_max, double b_fixed, int n, double* x_stable_upper, double* x_stable_lower, double* x_unstable, int* n_su, int* n_sl, int* n_unst);

/* Divergence analysis */
DivergenceTest* zeeman_divergence_test(double a_param, double b_start, double b_end, int n_trajectories);
void            zeeman_divergence_test_free(DivergenceTest* dt);
void            zeeman_divergence_test_print(const DivergenceTest* dt);

/* Inaccessibility & classification */
bool   zeeman_is_inaccessible(double a, double b, double x);
void   zeeman_map_to_cusp_normal(double c3, double c2, double c1, double c0, double* a_out, double* b_out, double* shift);

/* Potential surface */
double* zeeman_potential_surface(double a_min, double a_max, int na, double b_min, double b_max, int nb, double x);
double  zeeman_global_minimum(double a, double b);
double  zeeman_barrier_height(double a, double b);

/* Slow-fast cusp dynamics */
SlowFastCusp* slowfast_cusp_create(double epsilon, double a0, double b0, double alpha, double beta, double gamma);
void          slowfast_cusp_free(SlowFastCusp* sf);
void          slowfast_cusp_rhs(const SlowFastCusp* sf, double* dx, double* dy);
void          slowfast_cusp_step(SlowFastCusp* sf, double dt);
void          slowfast_cusp_simulate(SlowFastCusp* sf, double duration, double dt);
bool          slowfast_cusp_is_relaxation(const SlowFastCusp* sf);
double        slowfast_cusp_period(SlowFastCusp* sf, double duration, double dt);
void          slowfast_cusp_print(const SlowFastCusp* sf);

/* Zeeman's five qualities of catastrophe */
bool   zeeman_quality_bimodality(double a, double b);
bool   zeeman_quality_sudden_jump(const double* a, const double* x, int n, double* jump_a, double* jump_magnitude);
double zeeman_quality_hysteresis(double a_min, double a_max, double b, int steps);
double zeeman_quality_divergence(double a, double b_start, double b_end, int n);
bool   zeeman_quality_inaccessibility(double a, double b);
void   zeeman_five_qualities_report(double a, double b);

#endif /* ZEEMAN_DYNAMICS_H */