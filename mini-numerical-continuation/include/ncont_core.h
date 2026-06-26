#ifndef NCONT_CORE_H
#define NCONT_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * ncont_core.h ? Numerical Continuation Core Types
 *
 * Numerical continuation traces solution branches of
 *   F(x, lambda) = 0,  x in R^n, lambda in R
 *
 * as the parameter lambda varies. This is the computational
 * backbone of bifurcation analysis.
 *
 * Methods:
 *   - Natural parameter continuation
 *   - Pseudo-arclength continuation (Keller, 1977)
 *   - Moore-Penrose continuation
 *
 * Key concepts:
 *   - Predictor: extrapolate next point on the branch
 *   - Corrector: converge to the solution manifold
 *   - Step size control: adapt based on convergence
 *   - Bifurcation detection: test functions along the branch
 *
 * References:
 *   Keller (1977) Numerical Solution of Bifurcation Problems
 *   Allgower & Georg (2003) Numerical Continuation Methods
 *   Seydel (2010) Practical Bifurcation and Stability Analysis
 *   Kuznetsov (2004) Elements of Applied Bifurcation Theory
 * ============================================================== */

#define NCONT_MAX_DIM 20
#define NCONT_MAX_STEPS 10000

/* System function: F(x, lambda, params) -> residual */
typedef void (*NCONT_System)(const double* x, int n, double lambda,
                              void* params, double* residual);

/* Jacobian: dF/dx at (x, lambda) stored row-major */
typedef void (*NCONT_Jacobian)(const double* x, int n, double lambda,
                                void* params, double* jacobian);

/* Solution point on a branch */
typedef struct {
    double x[NCONT_MAX_DIM];
    int n;
    double lambda;
    double* tangent;         /* Unit tangent to the branch */
    double arclength;        /* Arclength parameter s */
    int step_index;
    bool converged;
    int newton_iterations;
    double residual_norm;
} NCONT_Point;

/* Continuation configuration */
typedef struct {
    double ds;               /* Initial step size */
    double ds_min;           /* Minimum step size */
    double ds_max;           /* Maximum step size */
    double newton_tol;       /* Newton convergence tolerance */
    int newton_max_iter;     /* Maximum Newton iterations */
    double arclength_max;    /* Maximum arclength to compute */
    int max_points;          /* Maximum number of points */
    bool detect_bifurcations;
    bool adaptive_stepsize;
    double step_adapt_factor; /* Factor for step size adaptation */
    int predictor_order;      /* 0=constant, 1=tangent, 2=secant */
} NCONT_Config;

/* Continuation state */
typedef struct {
    NCONT_Point* branch;     /* Array of solution points */
    int n_points;             /* Current number of points */
    int capacity;             /* Allocated capacity */
    NCONT_Config config;      /* Configuration */
    double* workspace1;       /* Workspace for linear algebra (n x n) */
    double* workspace2;       /* Additional workspace */
    double* rhs;              /* Right-hand side vector */
    int n;                    /* System dimension */
    bool finished;
    int exit_flag;            /* 0=success, 1=max_points, 2=arclength, 3=failure */
} NCONT_State;

/* Linear algebra utilities */
void ncont_solve_linear(double* A, double* b, int n, double* x);
bool ncont_solve_linear_safe(double* A, double* b, int n, double* x);
double ncont_vector_norm(const double* v, int n);
double ncont_vector_dot(const double* a, const double* b, int n);
void ncont_vector_copy(const double* src, double* dst, int n);
void ncont_vector_scale(double* v, int n, double s);
void ncont_vector_axpy(double* v, const double* w, int n, double a);

/* Solution point management */
NCONT_Point ncont_point_create(int n, const double* x, double lambda);
NCONT_Point ncont_point_copy(const NCONT_Point* p);
void ncont_point_set_tangent(NCONT_Point* p, const double* tangent, int n);

/* Continuation state lifecycle */
NCONT_State* ncont_state_create(int n, const double* x0, double lambda0,
                                 const NCONT_Config* config);
void ncont_state_free(NCONT_State* state);
int ncont_state_add_point(NCONT_State* state, const NCONT_Point* point);
void ncont_state_print(const NCONT_State* state);

/* Default configuration */
NCONT_Config ncont_config_default(void);
void ncont_config_print(const NCONT_Config* config);

#endif
