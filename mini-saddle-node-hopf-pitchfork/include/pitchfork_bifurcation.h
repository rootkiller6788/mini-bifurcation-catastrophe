#ifndef PITCHFORK_BIFURCATION_H
#define PITCHFORK_BIFURCATION_H
#include "gst_core.h"
#include <stdbool.h>
/* Pitchfork Bifurcation. Normal form: dx/dt = r*x - x^3 (super) or r*x + x^3 (sub).
 * Requires Z2 symmetry: f(-x,r) = -f(x,r).
 * Supercritical: 1 stable FP -> 3 FPs (x=0 unstable, x=+-sqrt(r) stable).
 * Subcritical: discontinuous with hysteresis.
 * Eigenvalue: simple zero crossing. Codimension 1. */
typedef struct { double r_critical, x_critical; double pf_coeff; bool is_supercritical, has_z2_sym; double* param_range; int n_params; double* sym_branch; double* asym_plus; double* asym_minus; int n_branch; } PitchforkBifurcation;
PitchforkBifurcation* pf_create(void);
void pf_free(PitchforkBifurcation* pb);
int pf_detect(ODEFunc f, double* x0, double* params, int n, int param_idx, double p_min, double p_max, double p_step, int max_iter, double tol, PitchforkBifurcation* result);
bool pf_check_z2_symmetry(ODEFunc f, int n);
double pf_normal_form_coeff(ODEFunc f, const double* xc, double* params, int n, int pidx);
int pf_compute_branches(PitchforkBifurcation* pb, ODEFunc f, double* x0, double* params, int n, int pidx, double p_min, double p_max, int n_pts, int max_iter, double tol);
bool pf_is_pitchfork(const EigenSpectrum* e, ODEFunc f, int n);
void pf_print(const PitchforkBifurcation* pb);
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
/* @section pitchfork_bifurcation_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
