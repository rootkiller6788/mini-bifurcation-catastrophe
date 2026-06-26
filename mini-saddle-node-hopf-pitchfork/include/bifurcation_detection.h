#ifndef BIFURCATION_DETECTION_H
#define BIFURCATION_DETECTION_H
#include "gst_core.h"
#include "saddle_node.h"
#include "hopf_bifurcation.h"
#include "pitchfork_bifurcation.h"
#include <stdbool.h>
/* General bifurcation detection: scans parameter ranges to find
 * all local bifurcations. Tracks fixed points and eigenvalues.
 * Classifies each bifurcation based on eigenvalue crossing patterns. */
typedef enum { BIF_SN=0, BIF_HOPF=1, BIF_PITCHFORK=2, BIF_TRANS=3, BIF_NONE=4 } BifurcationType;
typedef struct { BifurcationType type; double param; double* fp; int n; EigenSpectrum* eigs; } BifurcationEvent;
typedef struct { BifurcationEvent* events; int n, cap; double* param_scan; int n_scan; SaddleNodeBifurcation* sn; HopfBifurcation* hb; PitchforkBifurcation* pb; } BifurcationScanner;
BifurcationScanner* bscan_create(void);
void bscan_free(BifurcationScanner* bs);
int bscan_add_event(BifurcationScanner* bs, BifurcationType t, double p, double* fp, int n);
int bscan_detect_all(BifurcationScanner* bs, ODEFunc f, double* x0, double* params, int n, int param_idx, double p_min, double p_max, int n_steps, int max_iter, double tol);
int bscan_count_type(const BifurcationScanner* bs, BifurcationType t);
BifurcationEvent* bscan_get_event(const BifurcationScanner* bs, int idx);
const char* bscan_type_name(BifurcationType t);
void bscan_print(const BifurcationScanner* bs);
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
/* @section bifurcation_detection_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
