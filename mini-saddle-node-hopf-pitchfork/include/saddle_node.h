#ifndef SADDLE_NODE_H
#define SADDLE_NODE_H
#include "gst_core.h"
#include <stdbool.h>
/* Saddle-Node (Fold) Bifurcation.
 * Normal form: dx/dt = r + x^2 (or r - x^2).
 * r<0: two FPs (one stable, one unstable). r=0: one semi-stable FP.
 * r>0: no FPs. Eigenvalue: simple zero, others Re<0.
 * Codimension 1. Unfolding: r. */
typedef struct { double r_critical, x_critical; double eig_real, eig_imag; bool is_fold, is_supercritical; double* param_range; int n_params; double* stable_branch; double* unstable_branch; int n_branch; } SaddleNodeBifurcation;
SaddleNodeBifurcation* sn_create(void);
void sn_free(SaddleNodeBifurcation* sn);
int sn_detect(ODEFunc f, double* x0, double* params, int n, int param_idx, double p_min, double p_max, double p_step, int max_iter, double tol, SaddleNodeBifurcation* result);
int sn_compute_branches(SaddleNodeBifurcation* sn, ODEFunc f, double* x0, double* params, int n, int param_idx, double p_min, double p_max, int n_pts, int max_iter, double tol);
double sn_fold_condition(const Matrix* jac);
bool sn_is_saddle_node(const EigenSpectrum* e);
double sn_normal_form_coeff(ODEFunc f, const double* xc, double* params, int n, int pidx);
void sn_print(const SaddleNodeBifurcation* sn);
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
/* @section saddle_node_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
