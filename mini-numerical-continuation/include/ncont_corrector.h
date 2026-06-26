#ifndef NCONT_CORRECTOR_H
#define NCONT_CORRECTOR_H
#include "ncont_core.h"
#include "ncont_predictor.h"

/* ==============================================================
 * ncont_corrector.h ? Corrector Methods
 *
 * Given a predictor guess, converge to the solution manifold.
 * Uses Newton's method with pseudo-arclength constraint:
 *
 *   F(x, lambda) = 0
 *   (x - x_pred) dot t_x + (lambda - lambda_pred) * t_lambda - ds = 0
 *
 * or Moore-Penrose: step orthogonal to tangent.
 *
 * References:
 *   Keller (1977)
 *   Deuflhard (2004) Newton Methods for Nonlinear Problems
 * ============================================================== */

/* Corrector result */
typedef struct {
    double x[NCONT_MAX_DIM];
    int n;
    double lambda;
    bool converged;
    int iterations;
    double residual_norm;
    double step_norm;
} NCONT_Corrector;

/* Newton corrector with pseudo-arclength constraint */
NCONT_Corrector ncont_correct_arclength(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda, double ds,
    double tol, int max_iter);

/* Newton corrector with natural parameter (fix lambda) */
NCONT_Corrector ncont_correct_natural(NCONT_System F, NCONT_Jacobian J,
    const double* x0, int n, double lambda, void* params,
    double tol, int max_iter);

/* Newton corrector with Moore-Penrose (minimum norm) */
NCONT_Corrector ncont_correct_moore_penrose(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda,
    double tol, int max_iter);

/* Damped Newton with line search */
NCONT_Corrector ncont_correct_damped(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Predictor* pred, void* params,
    const double* tangent, double tangent_lambda, double ds,
    double tol, int max_iter, double damping_factor);

/* Check if corrector converged */
bool ncont_corrector_converged(const NCONT_Corrector* corr);

/* Estimate condition number at solution */
double ncont_estimate_condition(NCONT_Jacobian J, const double* x, int n,
                                 double lambda, void* params);

#endif

























