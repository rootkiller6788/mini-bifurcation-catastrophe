#ifndef BCD_WASHOUT_H
#define BCD_WASHOUT_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int bcd_washout_op0(int n, const double *x, double *out);
int bcd_washout_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int bcd_washout_op1(int n, const double *x, double *out);
int bcd_washout_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int bcd_washout_op2(int n, const double *x, double *out);
int bcd_washout_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int bcd_washout_op3(int n, const double *x, double *out);
int bcd_washout_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int bcd_washout_op4(int n, const double *x, double *out);
int bcd_washout_batch4(int n, int m, const double *x, double *out);
/* hyperbolic tangent element-wise transform */
int bcd_washout_op5(int n, const double *x, double *out);
int bcd_washout_batch5(int n, int m, const double *x, double *out);
double bcd_washout_mean(int n, const double *x);
double bcd_washout_variance(int n, const double *x);
int bcd_washout_minmax(int n, const double *x, double *out);
double bcd_washout_sum(int n, const double *x);
double bcd_washout_norm_l2(int n, const double *x);
double bcd_washout_norm_l1(int n, const double *x);
double bcd_washout_norm_linf(int n, const double *x);

#endif