#ifndef BCD_STABILITY_H
#define BCD_STABILITY_H
#include "bcd_core.h"

typedef struct { double* A; double* B; int n; double tau; double* eigenvalues; int n_eig; BCDStability stability; } BCDDelaySystem;
typedef struct { double* K; int n; double tau; double* poles; int n_poles; int is_stabilized; } BCDPolePlacement;
typedef struct { double* taus; double* max_re; int n; int switch_count; double* switch_taus; } BCDStabilityChart;

BCDDelaySystem* bcd_delay_system_create(const double* A, const double* B, int n, double tau);
void bcd_delay_system_free(BCDDelaySystem* ds);
BCDPolePlacement* bcd_pole_placement_design(const BCDDelaySystem* ds, const double* poles, int np);
void bcd_pole_placement_free(BCDPolePlacement* pp);
BCDDelaySystem* bcd_linearize_dde(const BCDDDE* dde, BCDDynamics f, void* p, const double* eq);
void bcd_delay_system_eigenvalues(BCDDelaySystem* ds, int n_terms);
BCDDelaySystem* bcd_stability_chart(const BCDDelaySystem* ds, double tmin, double tmax, int n_pts);
void bcd_stability_chart_free(BCDDelaySystem* sc);
int bcd_num_unstable_roots(const BCDDelaySystem* ds);
double bcd_stability_radius(const BCDDelaySystem* ds);
void bcd_stability_lobes(double tau, const double* A, const double* B, int n, double* lobe_taus, int max_lobes);
int bcd_stability_lobe_count(double tau_max, const double* A, const double* B, int n);
double bcd_stability_delay_margin(const BCDDelaySystem* ds);
int bcd_routh_hurwitz_delay(const double* poly, int deg, double tau, int* stable);
double bcd_lambert_w_principal(double x);
double bcd_min_stable_delay(const double* A, const double* B, int n, double search_max);
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
