#ifndef NCONT_ARCLENGTH_H
#define NCONT_ARCLENGTH_H
#include "ncont_core.h"
#include "ncont_predictor.h"
#include "ncont_corrector.h"

/* ==============================================================
 * ncont_arclength.h ? Pseudo-Arclength Continuation
 *
 * The main continuation loop:
 *   while (s < s_max && n_points < max_points):
 *       1. Predict next point using tangent/secant
 *       2. Correct to solution manifold
 *       3. Adapt step size
 *       4. Store point
 *       5. Check for bifurcations
 *
 * References:
 *   Keller (1977)
 *   Allgower & Georg (2003)
 * ============================================================== */

/* Continuation result summary */
typedef struct {
    int total_points;
    int accepted_points;
    int rejected_steps;
    double arclength_covered;
    double min_stepsize_used;
    double max_stepsize_used;
    double avg_newton_iterations;
    int n_bifurcations_detected;
    bool completed;
    int exit_flag;
} NCONT_Result;

/* Run full continuation */
NCONT_Result ncont_continue(NCONT_System F, NCONT_Jacobian J,
    const double* x0, int n, double lambda0,
    const NCONT_Config* config, void* params,
    NCONT_State** state_out);

/* Single continuation step */
bool ncont_step(NCONT_System F, NCONT_Jacobian J,
    NCONT_State* state, void* params);

/* Step size adaptation */
double ncont_adapt_stepsize(const NCONT_State* state, double current_ds,
    int newton_iterations, bool accepted);

/* Detect turning point (fold) */
bool ncont_detect_turning_point(const NCONT_State* state);

/* Print continuation progress */
void ncont_result_print(const NCONT_Result* result);

#endif

































