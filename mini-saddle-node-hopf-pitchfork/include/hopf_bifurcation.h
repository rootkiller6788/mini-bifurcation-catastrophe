#ifndef HOPF_BIFURCATION_H
#define HOPF_BIFURCATION_H
#include "gst_core.h"
#include <stdbool.h>
/* Hopf Bifurcation. Normal form: dz/dt = (r+i*omega)*z - (alpha+i*beta)*|z|^2*z.
 * r<0: stable FP. r=0: complex conj eigenvalues cross imag axis.
 * r>0: limit cycle (supercritical l1<0) or destroyed (subcritical l1>0).
 * First Lyapunov coefficient l1 determines type. Codimension 1. */
typedef struct { double r_critical, omega_critical; double first_lyapunov; double* lc_amplitude; int n_amp; double* param_range; int n_params; bool is_supercritical; EigenSpectrum* crit_eigs; double* jac_at_crit; int jac_size; } HopfBifurcation;
HopfBifurcation* hopf_create(void);
void hopf_free(HopfBifurcation* hb);
int hopf_detect(ODEFunc f, double* x0, double* params, int n, int param_idx, double p_min, double p_max, double p_step, int max_iter, double tol, HopfBifurcation* result);
bool hopf_is_hopf(const EigenSpectrum* e, double tol);
double hopf_first_lyapunov(ODEFunc f, const double* xc, double* params, int n, int pidx, double omega, double eps);
double hopf_lc_amplitude(const HopfBifurcation* hb, double r);
int hopf_detect_lc(ODEFunc f, double* x0, double* params, int n, double r, double dt, int n_periods, int n_trans, double* amp_out, double* period_out);
void hopf_print(const HopfBifurcation* hb);
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
/* @section hopf_bifurcation_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
