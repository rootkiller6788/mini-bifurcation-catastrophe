#ifndef NCONT_BIFURCATION_H
#define NCONT_BIFURCATION_H
#include "ncont_core.h"

/* ==============================================================
 * ncont_bifurcation.h ? Bifurcation Detection
 *
 * Test functions detect bifurcations along the branch:
 *   - Limit point (fold): d(lambda)/ds = 0
 *   - Branch point (pitchfork/transcritical): determinant sign change
 *   - Hopf point: pair of eigenvalues cross imaginary axis
 *
 * References:
 *   Seydel (2010) Ch.5
 *   Kuznetsov (2004) Ch.10
 * ============================================================== */

typedef enum {
    NCONT_BIF_NONE = 0,
    NCONT_BIF_LIMIT_POINT = 1,
    NCONT_BIF_BRANCH_POINT = 2,
    NCONT_BIF_HOPF = 3,
    NCONT_BIF_PERIOD_DOUBLING = 4
} NCONT_BifurcationType;

typedef struct {
    NCONT_BifurcationType type;
    double lambda;
    double x[NCONT_MAX_DIM];
    int n;
    int point_index;
    double test_value;
    double determinant_before;
    double determinant_after;
} NCONT_Bifurcation;

/* Test functions */
double ncont_test_lambda_turning(const NCONT_State* state, int idx);
double ncont_test_determinant(NCONT_Jacobian J, const NCONT_Point* point,
                               int n, double lambda, void* params);
double ncont_test_hopf(NCONT_Jacobian J, const NCONT_Point* point,
                        int n, double lambda, void* params);

/* Detect bifurcation between two consecutive points */
NCONT_Bifurcation ncont_detect_bifurcation(NCONT_Jacobian J,
    const NCONT_Point* p1, const NCONT_Point* p2,
    int n, void* params);

/* Scan branch for all bifurcations */
int ncont_scan_bifurcations(const NCONT_State* state, NCONT_Jacobian J,
    int n, void* params, NCONT_Bifurcation* bif_list, int max_bif);
void ncont_bifurcation_print(const NCONT_Bifurcation* bif);

/* Classify bifurcation */
const char* ncont_bifurcation_type_name(NCONT_BifurcationType type);
bool ncont_is_bifurcation_detected(const NCONT_Bifurcation* bif);

#endif

























