#ifndef TEC_UMBILICS_H
#define TEC_UMBILICS_H
#include "tec_core.h"

/* ==============================================================
 * tec_umbilics.h - Umbilic Catastrophes (D+4, D-4, D5)
 *
 * Two state variables (x, y):
 *   Hyperbolic (D+4): V = x^3 + y^3 + a*x*y + b*x + c*y
 *   Elliptic (D-4):   V = x^3 - 3*x*y^2 + a*(x^2+y^2) + b*x + c*y
 *   Parabolic (D5):   V = x^2*y + y^4 + a*x^2 + b*y^2 + c*x + d*y
 *
 * Gradient system: dx/dt = -dV/dx, dy/dt = -dV/dy
 * Equilibria: grad V = 0
 * Stability: eigenvalues of Hessian matrix
 * ============================================================== */

/* --- Hyperbolic Umbilic --- */
int tec_hyperbolic_equilibria(double a, double b, double c, double* xy_out, int max_pairs);
bool tec_hyperbolic_is_in_bifurcation_set(double a, double b, double c);
void tec_hyperbolic_phase_portrait(double a, double b, double c, double x0, double y0, double dt, int steps, double** traj_x, double** traj_y, int* n);

/* --- Elliptic Umbilic --- */
int tec_elliptic_equilibria(double a, double b, double c, double* xy_out, int max_pairs);
bool tec_elliptic_is_in_bifurcation_set(double a, double b, double c);
void tec_elliptic_phase_portrait(double a, double b, double c, double x0, double y0, double dt, int steps, double** traj_x, double** traj_y, int* n);

/* --- Parabolic Umbilic --- */
int tec_parabolic_equilibria(double a, double b, double c, double d, double* xy_out, int max_pairs);
bool tec_parabolic_is_in_bifurcation_set(double a, double b, double c, double d);
void tec_parabolic_phase_portrait(double a, double b, double c, double d, double x0, double y0, double dt, int steps, double** traj_x, double** traj_y, int* n);

/* --- Common Umbilic Utilities --- */
double tec_umbilic_potential_surface(TEC_CatastropheType type, const double* params, double x, double y);
void tec_umbilic_gradient_dynamics(TEC_CatastropheType type, const double* params, double* x, double* y, double dt);

#endif
