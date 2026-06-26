#ifndef BCD_CORE_H
#define BCD_CORE_H
#include <stddef.h>
#include <stdbool.h>

/* bcd_core.h - Bifurcation Control with Delay: Core Types, DDE Solver, Lambert W
 * dx/dt = f(x(t), x(t-tau), mu) where tau is the delay.
 * References: Stepan (1989) Retarded Dynamical Systems; Pyragas (1992) DFC */

#define BCD_MAX_DELAY    100.0
#define BCD_MAX_HISTORY  100000
#define BCD_MAX_DIM      32
#define BCD_PI           3.141592653589793
#define BCD_E            2.718281828459045
#define BCD_SQRT2        1.414213562373095
#define BCD_EPS          1e-10
#define BCD_DEFAULT_DT   0.01
#define BCD_DEFAULT_INTEGRATOR 0
#define BCD_GOLDEN       1.618033988749895
#define BCD_SQRT3        1.732050807568877

typedef enum { BCD_BIF_NONE=0, BCD_BIF_HOPF, BCD_BIF_SADDLE_NODE, BCD_BIF_PITCHFORK, BCD_BIF_PERIOD_DOUBLING, BCD_BIF_TORUS } BCDBifType;
typedef enum { BCD_STABLE=0, BCD_UNSTABLE, BCD_MARGINAL, BCD_UNDEFINED } BCDStability;
typedef enum { BCD_DDE_CONSTANT=0, BCD_DDE_DISTRIBUTED, BCD_DDE_STATE_DEPENDENT, BCD_DDE_NEUTRAL } BCDDDEType;

typedef void (*BCDDynamics)(int n, const double* x, const double* x_delayed, double t, double* dx, void* params);

typedef struct { int dim; double tau; double* history; int hist_len, hist_cap; double* buffer; int buf_pos, buf_cap; double t_current; BCDDDEType type; double* weights; int n_weights; void* extra; } BCDDDE;
typedef struct { double* x; double t; int step; } BCDDDEState;
typedef struct { double* eigenvalues; int n, has_hopf; double critical_tau, critical_freq; } BCDHopfResult;
typedef struct { double gain, delay; int dim; double* buffer; int buf_len, buf_pos; } BCDDFControl;
typedef struct { double omega, alpha; double* states; int n; } BCDWashoutFilter;
typedef struct { double* amps; int n_osc; double coupling, delay; int is_dead; } BCDAmplitudeDeath;

BCDDDE* bcd_dde_create(int dim, double tau, int hist_steps);
void bcd_dde_free(BCDDDE* dde);
void bcd_dde_set_history(BCDDDE* dde, const double* x0, int dim);
int bcd_dde_step(BCDDDE* dde, BCDDynamics f, void* params, double dt);
double bcd_dde_get_delayed(const BCDDDE* dde, int dim_idx);
void bcd_dde_simulate(BCDDDE* dde, BCDDynamics f, void* p, double T, double dt, double* out, int max_steps);
void bcd_dde_reset(BCDDDE* dde);

BCDDFControl* bcd_control_create(double gain, double delay, int dim);
void bcd_control_free(BCDDFControl* ctrl);
int bcd_control_apply(BCDDFControl* ctrl, const double* x, double t, double* u);
void bcd_control_update(BCDDFControl* ctrl, const double* x);

BCDWashoutFilter* bcd_washout_create(double omega, double alpha, int n);
void bcd_washout_free(BCDWashoutFilter* wf);
int bcd_washout_step(BCDWashoutFilter* wf, const double* x, double dt, double* y);

double bcd_lambert_w(double x, int branch);
int bcd_char_roots(double tau, const double* A, const double* B, int n, double* re, double* im, int max_roots);
int bcd_stability_switches(double tmin, double tmax, int n, const double* A, const double* B, int dim, double* sw, int max_sw);
double bcd_critical_delay(const double* A, const double* B, int n);

void bcd_dde_jacobian(BCDDDE* dde, BCDDynamics f, void* p, const double* x0, double eps, double* J);
void bcd_dde_delayed_jacobian(BCDDDE* dde, BCDDynamics f, void* p, const double* x0, double eps, double* Jd);

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
