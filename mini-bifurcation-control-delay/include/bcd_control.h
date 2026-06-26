#ifndef BCD_CONTROL_H
#define BCD_CONTROL_H
#include "bcd_core.h"

typedef enum { BCD_CTRL_PYRAGAS=0, BCD_CTRL_EXTENDED_DFC, BCD_CTRL_ADAPTIVE, BCD_CTRL_OPTIMAL, BCD_CTRL_PREDICTIVE } BCDControlType;
typedef struct { BCDControlType type; double gain, delay; double* target; int dim; double adapt_rate; } BCDControlConfig;
typedef struct { double* u_history; int n; double energy, settling_time, steady_error; int success; } BCDControlResult;

BCDControlConfig bcd_control_config_default(void);
BCDControlResult* bcd_control_design(const BCDDDE* dde, BCDDynamics f, void* p, const BCDControlConfig* cfg);
void bcd_control_result_free(BCDControlResult* cr);
double bcd_pyragas_optimal_gain(const BCDDDE* dde, BCDDynamics f, void* p, double omega);
double bcd_extended_dfc_gain(const BCDDDE* dde, BCDDynamics f, void* p, int n_terms);
int bcd_adaptive_dfc_step(BCDDFControl* ctrl, const double* x, double error, double rate);
double bcd_control_energy_efficiency(const BCDControlResult* cr);
int bcd_control_is_stabilized(const BCDControlResult* cr, double tol);
double bcd_nyquist_criterion_delay(const double* A, const double* B, int n, double tau);
void bcd_control_bode(const BCDDFControl* ctrl, double wmin, double wmax, int n, double* mag, double* phase);
double bcd_control_gain_margin(const BCDDFControl* ctrl, double omega);
double bcd_control_phase_margin(const BCDDFControl* ctrl, double omega);
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
