#ifndef NF_POINCARE_H
#define NF_POINCARE_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int nf_poincare_op0(int n, const double *x, double *out);
int nf_poincare_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int nf_poincare_op1(int n, const double *x, double *out);
int nf_poincare_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int nf_poincare_op2(int n, const double *x, double *out);
int nf_poincare_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int nf_poincare_op3(int n, const double *x, double *out);
int nf_poincare_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int nf_poincare_op4(int n, const double *x, double *out);
int nf_poincare_batch4(int n, int m, const double *x, double *out);
double nf_poincare_mean(int n, const double *x);
double nf_poincare_variance(int n, const double *x);
int nf_poincare_minmax(int n, const double *x, double *out);
double nf_poincare_sum(int n, const double *x);
double nf_poincare_norm_l2(int n, const double *x);
double nf_poincare_norm_l1(int n, const double *x);
double nf_poincare_norm_linf(int n, const double *x);

#endif