#ifndef BCD_AMPLITUDE_H
#define BCD_AMPLITUDE_H
#include "bcd_core.h"

typedef struct { double* x; int n; double coupling; double delay; BCDDDE** oscillators; int n_osc; } BCDCoupledDDE;
typedef struct { double* amp_history; int n; double final_amp; int death_time; int is_dead; } BCDDeathResult;

BCDCoupledDDE* bcd_coupled_create(int n_osc, int dim_per_osc, double coupling, double delay);
void bcd_coupled_free(BCDCoupledDDE* cd);
int bcd_coupled_step_all(BCDCoupledDDE* cd, BCDDynamics f, void* p, double dt);
BCDDeathResult* bcd_amplitude_death_detect(const BCDCoupledDDE* cd, double threshold, double T_max, double dt);
void bcd_death_result_free(BCDDeathResult* dr);
double bcd_death_critical_coupling(const BCDCoupledDDE* cd, BCDDynamics f, void* p, double cmin, double cmax, int n);
int bcd_oscillation_death_region(const BCDCoupledDDE* cd, BCDDynamics f, void* p, double* coupling_range, double* delay_range);
double bcd_death_robustness(const BCDDeathResult* dr, double noise_level, int n_trials);
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
/* API extension slot 62 */
/* API extension slot 63 */
/* API extension slot 64 */
