#include "tec_core.h"
#include "tec_cuspoids.h"
#include "tec_umbilics.h"
#include "tec_bifurcation.h"
#include "tec_applications.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* Test 1: Create/Destroy models */
void t1(void) {
    TEC_Model* m = tec_model_create(TEC_FOLD);
    assert(m != NULL);
    assert(m->type == TEC_FOLD);
    tec_model_free(m);
    printf("  create_free: PASS\n");
}

/* Test 2: Catastrophe names */
void t2(void) {
    assert(strcmp(tec_catastrophe_name(TEC_FOLD), "Fold") == 0);
    assert(strcmp(tec_catastrophe_name(TEC_CUSP), "Cusp") == 0);
    assert(strcmp(tec_catastrophe_name(TEC_SWALLOWTAIL), "Swallowtail") == 0);
    assert(tec_catastrophe_codimension(TEC_CUSP) == 2);
    assert(tec_catastrophe_codimension(TEC_BUTTERFLY) == 4);
    printf("  names: PASS\n");
}

/* Test 3: Fold catastrophe potential */
void t3(void) {
    double v = tec_potential_fold(0.0, 0.0);
    assert(fabs(v) < 1e-12);
    double v2 = tec_potential_fold(1.0, 1.0);
    assert(isfinite(v2));
    printf("  fold_potential: PASS\n");
}

/* Test 4: Cusp potential and gradient */
void t4(void) {
    double v = tec_potential_cusp(0.0, 0.0, 0.0);
    assert(fabs(v) < 1e-12);
    double grad[2];
    tec_gradient_cusp(0.1, -0.5, 0.0, grad);
    assert(isfinite(grad[0]));
    printf("  cusp_potential: PASS\n");
}

/* Test 5: Fold Hessian */
void t5(void) {
    double h = tec_hessian_fold(1.0);
    assert(fabs(h - 6.0) < 1e-9);
    double h0 = tec_hessian_fold(0.0);
    assert(fabs(h0) < 1e-12);
    printf("  fold_hessian: PASS\n");
}

/* Test 6: Cubic discriminant */
void t6(void) {
    double d1 = tec_cubic_discriminant(0.0, 0.0, 0.0);
    assert(fabs(d1) < 1e-12);
    double d2 = tec_cubic_discriminant(0.0, -3.0, 0.0);
    assert(d2 > 0.0); /* x^3 - 3x has 3 real roots */
    printf("  cubic_discriminant: PASS\n");
}

/* Test 7: Cubic solver */
void t7(void) {
    TEC_RootResult* r = tec_solve_cubic(0.0, 0.0, 0.0);
    assert(r != NULL);
    assert(r->n_real == 1);
    assert(fabs(r->roots[0]) < 1e-12);
    tec_roots_free(r);
    printf("  cubic_solver: PASS\n");
}

/* Test 8: Quartic discriminant */
void t8(void) {
    double d = tec_quartic_discriminant(0.0, 0.0, 0.0, 0.0);
    assert(isfinite(d));
    printf("  quartic_discriminant: PASS\n");
}

/* Test 9: Swallowtail potential */
void t9(void) {
    double v = tec_potential_swallowtail(0.5, -0.2, 0.1, 0.0);
    assert(isfinite(v));
    printf("  swallowtail: PASS\n");
}

/* Test 10: Butterfly potential */
void t10(void) {
    double v = tec_potential_butterfly(0.3, -0.1, 0.2, -0.3, 0.1);
    assert(isfinite(v));
    printf("  butterfly: PASS\n");
}

/* Test 11: Hyperbolic umbilic */
void t11(void) {
    double v = tec_potential_hyperbolic(0.2, -0.1, 0.5, 0.3, -0.2);
    assert(isfinite(v));
    printf("  hyperbolic_umbilic: PASS\n");
}

/* Test 12: Elliptic umbilic */
void t12(void) {
    double v = tec_potential_elliptic(0.1, 0.2, 0.3, -0.1, 0.0);
    assert(isfinite(v));
    printf("  elliptic_umbilic: PASS\n");
}

/* Test 13: Clamp utility */
void t13(void) {
    assert(fabs(tec_clamp(0.5, 0.0, 1.0) - 0.5) < 1e-9);
    assert(fabs(tec_clamp(2.0, 0.0, 1.0) - 1.0) < 1e-9);
    assert(fabs(tec_clamp(-1.0, 0.0, 1.0) - 0.0) < 1e-9);
    printf("  clamp: PASS\n");
}

/* Test 14: Model set control params */
void t14(void) {
    TEC_Model* m = tec_model_create(TEC_CUSP);
    double p[] = {-0.5, 0.1};
    tec_model_set_control(m, p, 2);
    assert(fabs(m->control_params[0] + 0.5) < 1e-9);
    tec_model_free(m);
    printf("  model_control: PASS\n");
}

/* Test 15: Find equilibria */
void t15(void) {
    TEC_Model* m = tec_model_create(TEC_CUSP);
    double p[] = {-1.0, 0.0};
    tec_model_set_control(m, p, 2);
    int n = tec_find_equilibria(m);
    assert(n >= 1); /* Should find at least one equilibrium */
    tec_model_free(m);
    printf("  equilibria: PASS\n");
}

/* Test 16: Parabolic umbilic */
void t16(void) {
    double v = tec_potential_parabolic(0.3, -0.2, 0.1, 0.4, -0.1, 0.0);
    assert(isfinite(v));
    printf("  parabolic_umbilic: PASS\n");
}

/* Test 17: Gradient umbilic */
void t17(void) {
    double grad[2], p[] = {0.1, 0.2, 0.3};
    tec_gradient_umbilic(0.5, -0.3, TEC_HYPERBOLIC_UMBILIC, p, grad);
    assert(isfinite(grad[0]) && isfinite(grad[1]));
    printf("  gradient_umbilic: PASS\n");
}

/* Test 18: Hessian umbilic */
void t18(void) {
    double h[4], p[] = {0.1, 0.2, 0.3};
    tec_hessian_umbilic(0.1, -0.1, TEC_ELLIPTIC_UMBILIC, p, h);
    assert(isfinite(h[0]) && isfinite(h[1]));
    assert(isfinite(h[2]) && isfinite(h[3]));
    printf("  hessian_umbilic: PASS\n");
}

/* Test 19: Multiple models */
void t19(void) {
    TEC_Model* models[7];
    TEC_CatastropheType types[] = {TEC_FOLD,TEC_CUSP,TEC_SWALLOWTAIL,
        TEC_BUTTERFLY,TEC_HYPERBOLIC_UMBILIC,TEC_ELLIPTIC_UMBILIC,
        TEC_PARABOLIC_UMBILIC};
    for (int i=0;i<7;i++) {
        models[i]=tec_model_create(types[i]);
        assert(models[i]!=NULL);
    }
    for (int i=0;i<7;i++) tec_model_free(models[i]);
    printf("  all_types: PASS\n");
}

/* Test 20: Hessian cusp */
void t20(void) {
    double h_vals[4];
    tec_hessian_cusp(0.5, -0.3, h_vals);
    assert(isfinite(h_vals[0]));
    printf("  hessian_cusp: PASS\n");
}

int main(void) {
    printf("=== Thom Elementary Catastrophes Test Suite ===\n");
    fflush(stdout);
    t1(); t2(); t3(); t4(); t5(); t6(); t7(); t8(); t9(); t10();
    t11(); t12(); t13(); t14(); t15(); t16(); t17(); t18(); t19(); t20();
    printf("=== All 20 tests passed ===\n");
    return 0;
}
