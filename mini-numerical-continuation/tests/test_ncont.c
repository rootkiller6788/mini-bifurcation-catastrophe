#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncont_core.h"
#include "ncont_predictor.h"
#include "ncont_corrector.h"
#include "ncont_arclength.h"
#include "ncont_bifurcation.h"
#include "ncont_branch.h"

#define EPS 1e-6

/* Simple test system: F(x,lambda) = x^2 - lambda */
void test_system(const double* x, int n, double lambda, void* params, double* residual) {
    (void)params;
    residual[0] = x[0] * x[0] - lambda;
}
void test_jacobian(const double* x, int n, double lambda, void* params, double* jac) {
    (void)params; (void)lambda;
    jac[0] = 2.0 * x[0];
}

void test_linear_algebra(void) {
    double A[4] = {2, 1, 1, 3}, b[2] = {5, 6}, x[2];
    ncont_solve_linear(A, b, 2, x);
    assert(fabs(x[0] - 1.8) < 0.1);
    assert(fabs(x[1] - 1.4) < 0.1);
}

void test_vector_ops(void) {
    double v[] = {3.0, 4.0};
    assert(fabs(ncont_vector_norm(v, 2) - 5.0) < EPS);
    assert(fabs(ncont_vector_dot(v, v, 2) - 25.0) < EPS);
}

void test_point_create(void) {
    double x[] = {1.0, 2.0};
    NCONT_Point p = ncont_point_create(2, x, 0.5);
    assert(p.n == 2);
    assert(fabs(p.x[0] - 1.0) < EPS);
    assert(fabs(p.lambda - 0.5) < EPS);
}

void test_state_create(void) {
    double x0[] = {1.0};
    NCONT_Config cfg = ncont_config_default();
    NCONT_State* s = ncont_state_create(1, x0, 2.0, &cfg);
    assert(s != NULL);
    assert(s->n_points == 1);
    assert(fabs(s->branch[0].x[0] - 1.0) < EPS);
    ncont_state_free(s);
}

void test_config_default(void) {
    NCONT_Config c = ncont_config_default();
    assert(c.ds > 0.0);
    assert(c.newton_max_iter > 0);
}

void test_tangent_predictor(void) {
    double x0[] = {1.0};
    NCONT_Config cfg = ncont_config_default(); cfg.predictor_order = 0;
    NCONT_State* s = ncont_state_create(1, x0, 2.0, &cfg);
    NCONT_Predictor pred = ncont_predict_tangent(s, 0, 0.1);
    assert(pred.n == 1);
    ncont_state_free(s);
}

void test_secant_predictor(void) {
    double x0[] = {1.0};
    NCONT_Config cfg = ncont_config_default();
    NCONT_State* s = ncont_state_create(1, x0, 2.0, &cfg);
    double x1[] = {1.1};
    NCONT_Point p1 = ncont_point_create(1, x1, 2.1);
    ncont_state_add_point(s, &p1);
    NCONT_Predictor pred = ncont_predict_secant(s, 0.1);
    assert(pred.n == 1);
    ncont_state_free(s);
}

void test_natural_corrector(void) {
    double x0[] = {1.5};
    NCONT_Corrector corr = ncont_correct_natural(test_system, test_jacobian, x0, 1, 2.25, NULL, 1e-6, 10);
    assert(corr.converged);
    assert(fabs(corr.x[0] - 1.5) < 0.01);
}

void test_arclength_corrector(void) {
    double x0[] = {1.5};
    NCONT_Config cfg = ncont_config_default();
    NCONT_State* s = ncont_state_create(1, x0, 2.25, &cfg);
    NCONT_Predictor pred = ncont_predict_tangent(s, 0, 0.1);
    double tangent[] = {1.0}; double tl = 0.5;
    NCONT_Corrector corr = ncont_correct_arclength(test_system, test_jacobian, &pred, NULL, tangent, tl, 0.1, 1e-6, 10);
    assert(corr.n == 1);
    ncont_state_free(s);
}

void test_continuation(void) {
    double x0[] = {1.0};
    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.05; cfg.arclength_max = 0.5; cfg.max_points = 20;
    NCONT_State* state = NULL;
    NCONT_Result result = ncont_continue(test_system, test_jacobian, x0, 1, 1.0, &cfg, NULL, &state);
    assert(state != NULL);
    assert(state->n_points >= 1);
    ncont_state_free(state);
}

void test_bifurcation_detection(void) {
    double x1[] = {0.1}, x2[] = {-0.1};
    NCONT_Point p1 = ncont_point_create(1, x1, 0.01);
    NCONT_Point p2 = ncont_point_create(1, x2, 0.01);
    p1.arclength = 0.0; p2.arclength = 0.1;
    NCONT_Bifurcation bif = ncont_detect_bifurcation(test_jacobian, &p1, &p2, 1, NULL);
    assert(bif.type == NCONT_BIF_LIMIT_POINT || bif.type == NCONT_BIF_NONE);
}

void test_branch_switch(void) {
    double xb[] = {0.0};
    NCONT_Point bp = ncont_point_create(1, xb, 0.0);
    NCONT_BranchSwitch bs = ncont_switch_branch_point(test_system, test_jacobian, &bp, 1, NULL, 0.1);
    assert(bs.success);
    assert(fabs(bs.x_new[0]) > 0.0);
}

void test_null_vector(void) {
    double x[] = {1.0, 0.0}, nv[2];
    bool ok = ncont_compute_null_vector(test_jacobian, x, 2, 1.0, NULL, nv);
    assert(ok);
}

void test_curve_handling(void) {
    double x0[] = {1.0, 0.5};
    NCONT_Config cfg = ncont_config_default();
    cfg.ds = 0.1; cfg.arclength_max = 0.3; cfg.max_points = 10;
    NCONT_State* s = ncont_state_create(2, x0, 1.0, &cfg);
    assert(s != NULL);
    ncont_state_free(s);
}

int main(void) {
    printf("=== Numerical Continuation Tests ===\n");
    test_linear_algebra();     printf("  linear_algebra PASSED\n");
    test_vector_ops();         printf("  vector_ops PASSED\n");
    test_point_create();       printf("  point_create PASSED\n");
    test_state_create();       printf("  state_create PASSED\n");
    test_config_default();     printf("  config PASSED\n");
    test_tangent_predictor();  printf("  tangent_predictor PASSED\n");
    test_secant_predictor();   printf("  secant_predictor PASSED\n");
    test_natural_corrector();  printf("  natural_corrector PASSED\n");
    test_arclength_corrector();printf("  arclength_corrector PASSED\n");
    test_continuation();       printf("  continuation PASSED\n");
    test_bifurcation_detection();printf("  bifurcation PASSED\n");
    test_branch_switch();      printf("  branch_switch PASSED\n");
    test_null_vector();        printf("  null_vector PASSED\n");
    test_curve_handling();     printf("  curve_handling PASSED\n");
    printf("=== All 14 tests passed ===\n");
    return 0;
}
