#ifndef TEC_CUSPOIDS_H
#define TEC_CUSPOIDS_H
#include "tec_core.h"

/* ==============================================================
 * tec_cuspoids.h - Cuspoid Catastrophes (A2, A3, A4, A5)
 *
 * One state variable, polynomial potentials:
 *   Fold (A2):       V = x^3 + a*x
 *   Cusp (A3):       V = x^4 + a*x^2 + b*x
 *   Swallowtail (A4): V = x^5 + a*x^3 + b*x^2 + c*x
 *   Butterfly (A5):  V = x^6 + a*x^4 + b*x^3 + c*x^2 + d*x
 *
 * Equilibrium: dV/dx = 0
 * Degeneracy: d^2V/dx^2 = 0
 *
 * Key surfaces:
 *   - Equilibrium manifold M: dV/dx = 0
 *   - Bifurcation set B: dV/dx = 0 AND d^2V/dx^2 = 0
 * ============================================================== */

/* --- Fold Catastrophe --- */
int tec_fold_equilibria(double a, double* x_out);
double tec_fold_bifurcation_point(void);
bool tec_fold_is_bistable(double a);
double tec_fold_jump_point(double a_start, double a_end, double x_start);

/* --- Cusp Catastrophe --- */
int tec_cusp_equilibria(double a, double b, double* x_out);
bool tec_cusp_is_in_bifurcation_set(double a, double b);
bool tec_cusp_is_bistable(double a, double b);
int tec_cusp_region_type(double a, double b);
double tec_cusp_maxwell_point(double a);
double tec_cusp_hysteresis_width(double a);

/* --- Swallowtail Catastrophe --- */
int tec_swallowtail_equilibria(double a, double b, double c, double* x_out);
bool tec_swallowtail_is_in_bifurcation_set(double a, double b, double c);

/* --- Butterfly Catastrophe --- */
int tec_butterfly_equilibria(double a, double b, double c, double d, double* x_out);
bool tec_butterfly_is_in_bifurcation_set(double a, double b, double c, double d);
int tec_butterfly_region_type(double a, double b, double c, double d);
double tec_butterfly_maxwell_set_symmetric(double b, double d);
double tec_butterfly_triple_point(double c);

double tec_cusp_energy_barrier(double a, double b);
double tec_cusp_bifurcation_curve_b(double a);
void tec_cusp_manifold_mesh(double am, double aM, double bm, double bM, int na, int nb, double** X, double** Y, double** Z);
int tec_cusp_region_analysis(double a, double b, double* xe, double* st, double* en);
double tec_fold_sensitivity(double a, double da);
double tec_cusp_sensitivity(double a, double b, double da, double db);
bool tec_butterfly_has_pocket(double a, double b, double c, double d);
double tec_butterfly_pocket_depth(double a, double b, double c, double d);
double tec_butterfly_pocket_boundary(double a, double b, double d, double* cr);
double tec_cusp_imperfection_sensitivity(double a, double b, double eps);

#endif
