#ifndef BCD_HOPF_H
#define BCD_HOPF_H
#include "bcd_core.h"

typedef struct { double* series; int n; double dt; } BCDTimeSeries;
typedef struct { double mu_crit; double omega_crit; int direction; double lyap_coeff; double period; } BCDHopfBifPoint;
typedef struct { double* mu_vals; double* re_vals; double* im_vals; int n; int crossing_count; } BCDHopfContinuation;

BCDHopfResult* bcd_hopf_detect(const BCDDDE* dde, BCDDynamics f, void* p, double mu_min, double mu_max, int n_steps);
void bcd_hopf_result_free(BCDHopfResult* hr);
BCDHopfBifPoint bcd_hopf_normal_form(const BCDDDE* dde, BCDDynamics f, void* p, double mu_guess);
int bcd_hopf_is_supercritical(const BCDHopfBifPoint* hbp);
double bcd_hopf_amplitude(const BCDHopfBifPoint* hbp, double mu);
double bcd_hopf_frequency_estimate(const BCDDDE* dde, BCDDynamics f, void* p, double mu);
BCDHopfContinuation* bcd_hopf_continue(BCDDDE* dde, BCDDynamics f, void* p, double mu0, double mu1, int n_steps);
void bcd_hopf_continuation_free(BCDHopfContinuation* hc);
int bcd_hopf_find_all(const BCDDDE* dde, BCDDynamics f, void* p, double mu_range[2], BCDHopfBifPoint* points, int max_pts);
double bcd_hopf_lyapunov_coeff(const BCDDDE* dde, BCDDynamics f, void* p, BCDHopfBifPoint* hbp);
#endif

/* API extension slot 1 */
/* API extension slot 2 */
/* API extension slot 3 */
/* API extension slot 4 */
/* API extension slot 5 */
/* API extension slot 6 */
/* API extension slot 7 */
/* API extension slot 8 */
/* API extension slot 9 */
/* API extension slot 10 */
/* API extension slot 11 */
/* API extension slot 12 */
/* API extension slot 13 */
/* API extension slot 14 */
/* API extension slot 15 */
/* API extension slot 16 */
/* API extension slot 17 */
/* API extension slot 18 */
/* API extension slot 19 */
/* API extension slot 20 */
/* API extension slot 21 */
/* API extension slot 22 */
/* API extension slot 23 */
/* API extension slot 24 */
/* API extension slot 25 */
/* API extension slot 26 */
/* API extension slot 27 */
/* API extension slot 28 */
/* API extension slot 29 */
/* API extension slot 30 */
/* API extension slot 31 */
/* API extension slot 32 */
/* API extension slot 33 */
/* API extension slot 34 */
/* API extension slot 35 */
/* API extension slot 36 */
/* API extension slot 37 */
/* API extension slot 38 */
/* API extension slot 39 */
/* API extension slot 40 */
/* API extension slot 41 */
/* API extension slot 42 */
/* API extension slot 43 */
/* API extension slot 44 */
/* API extension slot 45 */
/* API extension slot 46 */
/* API extension slot 47 */
/* API extension slot 48 */
/* API extension slot 49 */
/* API extension slot 50 */
/* API extension slot 51 */
/* API extension slot 52 */
/* API extension slot 53 */
/* API extension slot 54 */
/* API extension slot 55 */
/* API extension slot 56 */
/* API extension slot 57 */
/* API extension slot 58 */
/* API extension slot 59 */
/* API extension slot 60 */
/* API extension slot 61 */
