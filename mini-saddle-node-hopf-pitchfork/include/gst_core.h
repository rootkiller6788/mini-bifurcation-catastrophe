#ifndef GST_CORE_H
#define GST_CORE_H
#include <stdbool.h>
#include <stddef.h>
typedef struct { double* data; int n; } Vector;
typedef struct { double* data; int rows; int cols; } Matrix;
typedef struct { double real, imag; bool is_stable; double magnitude; } Eigenvalue;
typedef struct { Eigenvalue* values; int n, cap; } EigenSpectrum;
typedef void (*ODEFunc)(const double* x, double* params, double* dx, int n);
typedef struct { double* point; int n; EigenSpectrum* eigenvalues; bool is_stable, is_hyperbolic; double* jac_data; int jac_rows; } FixedPoint;
typedef enum { NORM_SN=0, NORM_HOPF=1, NORM_PF_SUPER=2, NORM_PF_SUB=3, NORM_TRANS=4, NORM_NONE=5 } NormalFormType;
typedef struct { NormalFormType type; double* coeffs; int n_coeff; double bif_param; char* expr; } NormalForm;

Vector* vec_create(int n); void vec_free(Vector* v); double vec_norm(const Vector* v);
double vec_dot(const Vector* a, const Vector* b);
Vector* vec_add(const Vector* a, const Vector* b);
Vector* vec_scale(const Vector* v, double s);

Matrix* mat_create(int r, int c); void mat_free(Matrix* m);
Matrix* mat_identity(int n); Matrix* mat_copy(const Matrix* m);
Matrix* mat_mul(const Matrix* a, const Matrix* b);
Matrix* mat_add(const Matrix* A, const Matrix* B);
Matrix* mat_sub(const Matrix* A, const Matrix* B);
Matrix* mat_transpose(const Matrix* m);
Matrix* mat_inverse(const Matrix* m);
double mat_trace(const Matrix* m); double mat_det(const Matrix* m);
void mat_set(Matrix* m, int i, int j, double v);
double mat_get(const Matrix* m, int i, int j);

EigenSpectrum* eigen_create(void); void eigen_free(EigenSpectrum* e);
int eigen_add(EigenSpectrum* e, double real, double imag);
EigenSpectrum* eigen_compute_2x2(const Matrix* J);
EigenSpectrum* eigen_compute_jacobian(ODEFunc f, const double* x, double* params, int n, double eps);
bool eigen_all_stable(const EigenSpectrum* e);

FixedPoint* fp_create(int n); void fp_free(FixedPoint* fp);
int fp_newton(ODEFunc f, double* x0, double* params, int n, int max_iter, double tol, FixedPoint* result);

void rk4_step(ODEFunc f, double* x, double* params, int n, double dt);
void rk4_solve(ODEFunc f, double* x0, double* params, int n, double t0, double tf, double dt, double** traj, int* n_steps);

NormalForm* nf_create(NormalFormType type); void nf_free(NormalForm* nf);
const char* nf_type_name(NormalFormType type);
void sn_normal_form(const double* x, double* params, double* dx, int n);
void hopf_normal_form(const double* x, double* params, double* dx, int n);
void pitchfork_normal_form(const double* x, double* params, double* dx, int n);
#endif
/* Utility functions */
void validate_bifurcation_inputs(ODEFunc f, const double* x, double* params, int n);
double safe_divide(double a, double b);
int sign_change_detect(const double* values, int n);
double min_positive(const double* values, int n);

/* @note All functions in this header are thread-safe if given
 * separate state. Functions with non-const pointer arguments
 * modify data in place. See gst_core.h for base types. */
/* @section gst_core_extended Extended API
 * This header is part of the local bifurcation analysis framework.
 * All functions perform parameter validation before computation.
 * Memory ownership: create functions allocate, free functions release.
 * Thread safety: functions use no global state. */
