#ifndef BIF_CODIM2_H
#define BIF_CODIM2_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int bif_codim2_op0(int n, const double *x, double *out);
int bif_codim2_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int bif_codim2_op1(int n, const double *x, double *out);
int bif_codim2_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int bif_codim2_op2(int n, const double *x, double *out);
int bif_codim2_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int bif_codim2_op3(int n, const double *x, double *out);
int bif_codim2_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int bif_codim2_op4(int n, const double *x, double *out);
int bif_codim2_batch4(int n, int m, const double *x, double *out);
double bif_codim2_mean(int n, const double *x);
double bif_codim2_variance(int n, const double *x);
int bif_codim2_minmax(int n, const double *x, double *out);
double bif_codim2_sum(int n, const double *x);
double bif_codim2_norm_l2(int n, const double *x);
double bif_codim2_norm_l1(int n, const double *x);
double bif_codim2_norm_linf(int n, const double *x);

#endif