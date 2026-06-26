#ifndef HOMOCLINIC_SHILNIKOV_H
#define HOMOCLINIC_SHILNIKOV_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int homoclinic_shilnikov_op0(int n, const double *x, double *out);
int homoclinic_shilnikov_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int homoclinic_shilnikov_op1(int n, const double *x, double *out);
int homoclinic_shilnikov_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int homoclinic_shilnikov_op2(int n, const double *x, double *out);
int homoclinic_shilnikov_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int homoclinic_shilnikov_op3(int n, const double *x, double *out);
int homoclinic_shilnikov_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int homoclinic_shilnikov_op4(int n, const double *x, double *out);
int homoclinic_shilnikov_batch4(int n, int m, const double *x, double *out);
double homoclinic_shilnikov_mean(int n, const double *x);
double homoclinic_shilnikov_variance(int n, const double *x);
int homoclinic_shilnikov_minmax(int n, const double *x, double *out);
double homoclinic_shilnikov_sum(int n, const double *x);
double homoclinic_shilnikov_norm_l2(int n, const double *x);
double homoclinic_shilnikov_norm_l1(int n, const double *x);
double homoclinic_shilnikov_norm_linf(int n, const double *x);

#endif