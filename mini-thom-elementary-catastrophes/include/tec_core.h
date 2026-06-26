#ifndef TEC_CORE_H
#define TEC_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * tec_core.h - Thom Elementary Catastrophe Theory Core Types
 *
 * Rene Thom (1972) "Stabilite Structurelle et Morphogenese"
 * Classified all structurally stable discontinuous changes that
 * can occur in gradient systems with <= 4 parameters.
 *
 * Seven elementary catastrophes:
 *   Cuspoids (1 variable):
 *     Fold (A2):        V = x^3 + a*x             [2 params]
 *     Cusp (A3):        V = x^4 + a*x^2 + b*x     [3 params]
 *     Swallowtail (A4): V = x^5 + a*x^3 + b*x^2 + c*x [4 params]
 *     Butterfly (A5):   V = x^6 + a*x^4 + b*x^3 + c*x^2 + d*x [5 params]
 *   Umbilics (2 variables):
 *     Hyperbolic (D+4): V = x^3 + y^3 + a*x*y + b*x + c*y
 *     Elliptic (D-4):   V = x^3 - 3*x*y^2 + a*(x^2+y^2) + b*x + c*y
 *     Parabolic (D5):   V = x^2*y + y^4 + a*x^2 + b*y^2 + c*x + d*y
 *
 * Key concepts:
 *   - Potential function V(x; c): smooth function of state x and control c
 *   - Equilibrium manifold M: grad_x V = 0
 *   - Singularity set S: det(Hess V) = 0 (degenerate critical points)
 *   - Bifurcation set B: projection of S onto control space
 *   - Maxwell set: where two minima have equal depth
 *
 * References:
 *   Thom (1972) Stabilite Structurelle et Morphogenese, Benjamin
 *   Zeeman (1977) Catastrophe Theory: Selected Papers, Addison-Wesley
 *   Poston & Stewart (1978) Catastrophe Theory and its Applications
 *   Gilmore (1981) Catastrophe Theory for Scientists and Engineers
 * ============================================================== */

/* --- Catastrophe Type Enumeration --- */
typedef enum {
    TEC_FOLD = 0,
    TEC_CUSP = 1,
    TEC_SWALLOWTAIL = 2,
    TEC_BUTTERFLY = 3,
    TEC_HYPERBOLIC_UMBILIC = 4,
    TEC_ELLIPTIC_UMBILIC = 5,
    TEC_PARABOLIC_UMBILIC = 6
} TEC_CatastropheType;

/* --- Critical Point Classification --- */
typedef enum {
    TEC_MINIMUM = 0,       /* Stable equilibrium */
    TEC_MAXIMUM = 1,       /* Unstable equilibrium */
    TEC_SADDLE = 2,         /* Saddle point (umbilics) */
    TEC_DEGENERATE = 3,    /* Hessian zero (bifurcation point) */
    TEC_INFLECTION = 4     /* Inflection point */
} TEC_CriticalType;

/* --- Real Polynomial Root (for cubic/quartic) --- */
typedef struct {
    double real;
    double imag;
    bool is_real;
} TEC_ComplexRoot;

typedef struct {
    TEC_ComplexRoot* roots;
    int n_roots;
    int n_real;
} TEC_RootResult;

/* --- Critical Point --- */
typedef struct {
    double x;               /* State variable value */
    double y;               /* Second state (umbilics only) */
    double V;               /* Potential value at this point */
    TEC_CriticalType type;
    double hessian_det;     /* Determinant of Hessian matrix */
    double hessian_trace;   /* Trace of Hessian */
    bool is_degenerate;
} TEC_CriticalPoint;

/* --- Control Parameter Vector --- */
typedef struct {
    double* params;         /* Control parameters (a, b, c, d) */
    int n_params;           /* Number of control parameters */
} TEC_ControlParams;

/* --- Bifurcation Set Point --- */
typedef struct {
    double* control;        /* Control parameter values */
    int n_control;
    double* critical_x;     /* State at bifurcation */
    int n_critical;
    TEC_CatastropheType type;
} TEC_BifurcationPoint;

/* --- Maxwell Set Point (equal-depth minima) --- */
typedef struct {
    double* control;
    int n_control;
    double x1, x2;          /* Two equal-depth minima positions */
    double V_min;           /* Common depth */
} TEC_MaxwellPoint;

/* --- Catastrophe Model --- */
typedef struct {
    TEC_CatastropheType type;
    char* name;
    int n_state_vars;       /* 1 for cuspoids, 2 for umbilics */
    int n_control_params;   /* Codimension */
    TEC_ControlParams control;
    TEC_CriticalPoint* equilibria;
    int n_equilibria;
} TEC_Model;

/* ==============================================================
 * Core API: Model Management
 * ============================================================== */
TEC_Model* tec_model_create(TEC_CatastropheType type);
void tec_model_free(TEC_Model* m);
void tec_model_set_control(TEC_Model* m, const double* params, int n);
void tec_model_print(const TEC_Model* m);
const char* tec_catastrophe_name(TEC_CatastropheType type);
int tec_catastrophe_codimension(TEC_CatastropheType type);

/* ==============================================================
 * Core API: Potential Functions
 * ============================================================== */
double tec_potential_fold(double x, double a);
double tec_potential_cusp(double x, double a, double b);
double tec_potential_swallowtail(double x, double a, double b, double c);
double tec_potential_butterfly(double x, double a, double b, double c, double d);
double tec_potential_hyperbolic(double x, double y, double a, double b, double c);
double tec_potential_elliptic(double x, double y, double a, double b, double c);
double tec_potential_parabolic(double x, double y, double a, double b, double c, double d);

/* ==============================================================
 * Core API: Gradient and Hessian
 * ============================================================== */
double tec_gradient_fold(double x, double a);
void tec_gradient_cusp(double x, double a, double b, double* grad);
void tec_gradient_umbilic(double x, double y, TEC_CatastropheType type, const double* params, double* grad);
double tec_hessian_fold(double x);
void tec_hessian_cusp(double x, double a, double* h);
void tec_hessian_umbilic(double x, double y, TEC_CatastropheType type, const double* p, double h[4]);

/* ==============================================================
 * Core API: Root Finding (Cubic, Quartic)
 * ============================================================== */
TEC_RootResult* tec_solve_cubic(double a, double b, double c);
TEC_RootResult* tec_solve_quartic(double a, double b, double c, double d);
void tec_roots_free(TEC_RootResult* r);

/* ==============================================================
 * Core API: Critical Point Analysis
 * ============================================================== */
int tec_find_equilibria(TEC_Model* m);
void tec_classify_critical_point(TEC_CriticalPoint* cp, TEC_CatastropheType type, const double* p);
void tec_critical_free(TEC_CriticalPoint* cp, int n);

/* ==============================================================
 * Core Utility
 * ============================================================== */
double tec_cubic_discriminant(double a, double b, double c);
double tec_quartic_discriminant(double a, double b, double c, double d);
double tec_clamp(double x, double lo, double hi);

#endif
