#ifndef CENTER_MANIFOLD_H
#define CENTER_MANIFOLD_H
#include "gst_core.h"
#include <stdbool.h>
/* Center Manifold Reduction.
 * For systems at bifurcation (zero real-part eigenvalues),
 * reduces dynamics to center manifold where bifurcation occurs.
 * Splitting: stable + center + unstable eigenspaces.
 * Reduction: x_dot = f(x,y), y_dot = g(x,y) with y = h(x)
 *   on center manifold. Then reduced system: x_dot = f(x, h(x)).
 * Theorem: local dynamics of full system is topologically equivalent
 *   to dynamics on center manifold (Center Manifold Theorem). */
typedef struct { Matrix* stable_basis; Matrix* center_basis; Matrix* unstable_basis; int n_stable, n_center, n_unstable; int n_total; double* transformation; int trans_size; } CenterManifold;
CenterManifold* cm_create(int n);
void cm_free(CenterManifold* cm);
int cm_compute(ODEFunc f, const double* x_crit, double* params, int n, double eps, CenterManifold* cm);
int cm_reduce(ODEFunc f, const double* x_crit, double* params, int n, CenterManifold* cm, ODEFunc* reduced_f);
int cm_approximate(ODEFunc f, const double* x_crit, double* params, int n, int order, double* h_coeffs, int* n_coeffs);
void cm_print(const CenterManifold* cm);
#endif/* Utility functions for numerical analysis */
double parameter_sensitivity(ODEFunc f, const double* x, double* params, int n, int pidx, double eps);
double eigenvalue_sensitivity(const Matrix* J, int pidx, double eps);
void print_bifurcation_summary(void);
/* Utility functions */
void validate_bifurcation_inputs(ODEFunc f, const double* x, double* params, int n);
double safe_divide(double a, double b);
int sign_change_detect(const double* values, int n);
double min_positive(const double* values, int n);

/* @note All functions in this header are thread-safe if given
 * separate state. Functions with non-const pointer arguments
 * modify data in place. See gst_core.h for base types. */
/* @section center_manifold_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
