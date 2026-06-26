#ifndef NCONT_PREDICTOR_H
#define NCONT_PREDICTOR_H
#include "ncont_core.h"

/* ==============================================================
 * ncont_predictor.h ? Predictor Methods
 *
 * Given points on the branch, extrapolate to guess the next point.
 *
 * Tangent predictor:  x_pred = x_k + ds * t_k
 *   where t_k is the unit tangent vector at the current point.
 *
 * Secant predictor:   x_pred = x_k + ds * (x_k - x_{k-1}) / |x_k - x_{k-1}|
 *   uses the secant direction from previous two points.
 *
 * Polynomial predictor: uses higher-order extrapolation.
 *
 * References:
 *   Allgower & Georg (2003) Ch.4
 * ============================================================== */

/* Predictor result: contains both state and parameter */
typedef struct {
    double x_pred[NCONT_MAX_DIM];
    double lambda_pred;
    int n;
    double predictor_norm;
} NCONT_Predictor;

/* Tangent predictor */
NCONT_Predictor ncont_predict_tangent(const NCONT_State* state, int idx,
                                        double ds);
void ncont_compute_tangent(NCONT_System F, NCONT_Jacobian J,
    const double* x, int n, double lambda, void* params,
    double* tangent, double* tangent_lambda);

/* Secant predictor (from last two points) */
NCONT_Predictor ncont_predict_secant(const NCONT_State* state, double ds);

/* Polynomial predictor (degree 0, 1, or 2) */
NCONT_Predictor ncont_predict_polynomial(const NCONT_State* state,
    double ds, int degree);

/* Select best predictor based on configuration */
NCONT_Predictor ncont_predict(const NCONT_State* state, double ds);

/* Compute tangent via bordered matrix method */
bool ncont_compute_tangent_bordered(NCONT_System F, NCONT_Jacobian J,
    const double* x, int n, double lambda, void* params,
    double* tangent, double* tangent_lambda);

/* Estimate curvature of the branch for step size control */
double ncont_estimate_curvature(const NCONT_State* state, int idx);

#endif

























