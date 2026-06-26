#include "ncont_bifurcation.h"
#include "ncont_predictor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Test Function: d(lambda)/ds for turning points
 * ============================================================== */

double ncont_test_lambda_turning(const NCONT_State* state, int idx) {
    if (!state || idx < 1 || idx >= state->n_points) return 0.0;
    double ds = state->branch[idx].arclength - state->branch[idx - 1].arclength;
    if (fabs(ds) < 1e-15) return 0.0;
    return (state->branch[idx].lambda - state->branch[idx - 1].lambda) / ds;
}

/* ==============================================================
 * Test Function: Determinant of Jacobian
 * ============================================================== */

double ncont_test_determinant(NCONT_Jacobian J, const NCONT_Point* point,
                               int n, double lambda, void* params) {
    if (!J || !point || n <= 0) return 0.0;
    double* jac = (double*)malloc(n * n * sizeof(double));
    J(point->x, n, lambda, params, jac);
    double det = 1.0;
    if (n == 1) det = jac[0];
    else if (n == 2) det = jac[0] * jac[3] - jac[1] * jac[2];
    else if (n == 3) {
        double* d = jac;
        det = d[0]*(d[4]*d[8]-d[5]*d[7]) - d[1]*(d[3]*d[8]-d[5]*d[6]) + d[2]*(d[3]*d[7]-d[4]*d[6]);
    }
    free(jac);
    return det;
}

/* ==============================================================
 * Test Function: Hopf detection (trace of Jacobian sign change)
 * ============================================================== */

double ncont_test_hopf(NCONT_Jacobian J, const NCONT_Point* point,
                        int n, double lambda, void* params) {
    if (!J || !point || n < 2) return 0.0;
    double* jac = (double*)malloc(n * n * sizeof(double));
    J(point->x, n, lambda, params, jac);
    double trace = 0.0; int i;
    for (i = 0; i < n; i++) trace += jac[i * n + i];
    free(jac);
    return trace;
}

/* ==============================================================
 * Bifurcation Detection
 * ============================================================== */

NCONT_Bifurcation ncont_detect_bifurcation(NCONT_Jacobian J,
    const NCONT_Point* p1, const NCONT_Point* p2, int n, void* params) {
    NCONT_Bifurcation bif; memset(&bif, 0, sizeof(bif));
    if (!J || !p1 || !p2 || n <= 0) return bif;
    bif.n = n; bif.lambda = p2->lambda; memcpy(bif.x, p2->x, n * sizeof(double));

    double det1 = ncont_test_determinant(J, p1, n, p1->lambda, params);
    double det2 = ncont_test_determinant(J, p2, n, p2->lambda, params);
    bif.determinant_before = det1; bif.determinant_after = det2;

    double dlambda = p2->lambda - p1->lambda;
    if (fabs(dlambda) < 1e-10) return bif;
    double dl_ds = dlambda / fmax(p2->arclength - p1->arclength, 1e-10);

    if (det1 * det2 < 0 && fabs(dl_ds) < 0.1) {
        bif.type = NCONT_BIF_BRANCH_POINT; bif.test_value = det2;
    } else if (det1 * det2 < 0) {
        bif.type = NCONT_BIF_LIMIT_POINT; bif.test_value = dlambda;
    } else {
        double tr1 = ncont_test_hopf(J, p1, n, p1->lambda, params);
        double tr2 = ncont_test_hopf(J, p2, n, p2->lambda, params);
        if (tr1 * tr2 < 0) { bif.type = NCONT_BIF_HOPF; bif.test_value = tr2; }
    }
    return bif;
}

/* ==============================================================
 * Scan Branch for Bifurcations
 * ============================================================== */

int ncont_scan_bifurcations(const NCONT_State* state, NCONT_Jacobian J,
    int n, void* params, NCONT_Bifurcation* bif_list, int max_bif) {
    if (!state || !J || !bif_list || max_bif <= 0) return 0;
    int count = 0, i;
    for (i = 1; i < state->n_points && count < max_bif; i++) {
        NCONT_Bifurcation bif = ncont_detect_bifurcation(J, &state->branch[i - 1], &state->branch[i], n, params);
        if (bif.type != NCONT_BIF_NONE) {
            bif.point_index = i;
            bif_list[count++] = bif;
        }
    }
    return count;
}

void ncont_bifurcation_print(const NCONT_Bifurcation* bif) {
    if (!bif) return;
    const char* names[] = {"None","LimitPoint","BranchPoint","Hopf","PeriodDoubling"};
    printf("Bifurcation: %s at lambda=%.6f (point %d) test=%.6f det:%.4f->%.4f\n",
        names[bif->type < 5 ? bif->type : 0], bif->lambda, bif->point_index,
        bif->test_value, bif->determinant_before, bif->determinant_after);
}

const char* ncont_bifurcation_type_name(NCONT_BifurcationType type) {
    const char* names[] = {"None","LimitPoint","BranchPoint","Hopf","PeriodDoubling"};
    return names[type < 5 ? type : 0];
}

bool ncont_is_bifurcation_detected(const NCONT_Bifurcation* bif) {
    return bif && bif->type != NCONT_BIF_NONE;
}
bool ncont_eigenvalue_sign_change(NCONT_Jacobian J, const NCONT_Point* p1, const NCONT_Point* p2, int n, void* params) {
    if(!J||!p1||!p2||n<=0)return false;
    double d1=ncont_test_determinant(J,p1,n,p1->lambda,params);
    double d2=ncont_test_determinant(J,p2,n,p2->lambda,params);return d1*d2<0;
}
int ncont_bifurcation_classify_detailed(NCONT_Jacobian J,
    const NCONT_Point* p1, const NCONT_Point* p2, int n, void* params) {
    double d1=ncont_test_determinant(J,p1,n,p1->lambda,params);
    double d2=ncont_test_determinant(J,p2,n,p2->lambda,params);
    double dl=p2->lambda-p1->lambda;
    if(d1*d2<0&&fabs(dl)<1e-3)return NCONT_BIF_LIMIT_POINT;
    if(d1*d2<0)return NCONT_BIF_BRANCH_POINT;return NCONT_BIF_NONE;
}
void ncont_bifurcation_count_by_type(const NCONT_Bifurcation* bl, int nb, int c[5]) {
    int i;for(i=0;i<5;i++)c[i]=0;for(i=0;i<nb;i++)if(bl[i].type<5)c[bl[i].type]++;
}
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* */
/* bifurcation block 1 */
/* bifurcation block 2 */
/* bifurcation block 3 */
/* bifurcation block 4 */
/* bifurcation block 5 */
/* bifurcation block 6 */
/* bifurcation block 7 */
/* bifurcation block 8 */
/* bifurcation block 9 */
/* bifurcation block 10 */
/* bifurcation block 11 */
/* bifurcation block 12 */
/* bifurcation block 13 */
/* bifurcation block 14 */
/* bifurcation block 15 */
/* bifurcation block 16 */
/* bifurcation block 17 */
/* bifurcation block 18 */
/* bifurcation block 19 */
/* bifurcation block 20 */
/* bifurcation block 21 */
/* bifurcation block 22 */
/* bifurcation block 23 */
/* bifurcation block 24 */
/* bifurcation block 25 */
/* bifurcation block 26 */
/* bifurcation block 27 */
/* bifurcation block 28 */
/* bifurcation block 29 */
/* bifurcation block 30 */
/* bifurcation block 31 */
/* bifurcation block 32 */
/* bifurcation block 33 */
/* bifurcation block 34 */
/* bifurcation block 35 */
/* bifurcation block 36 */
/* bifurcation block 37 */
/* bifurcation block 38 */
/* bifurcation block 39 */
/* bifurcation block 40 */
/* bifurcation block 41 */
/* bifurcation block 42 */
/* bifurcation block 43 */
/* bifurcation block 44 */
/* bifurcation block 45 */
/* bifurcation block 46 */
/* bifurcation block 47 */
/* bifurcation block 48 */
/* bifurcation block 49 */
/* bifurcation block 50 */
/* bifurcation block 51 */
/* bifurcation block 52 */
/* bifurcation block 53 */
/* bifurcation block 54 */
/* bifurcation block 55 */
/* bifurcation block 56 */
/* bifurcation block 57 */
/* bifurcation block 58 */
/* bifurcation block 59 */
/* bifurcation block 60 */
/* bifurcation block 61 */
/* bifurcation block 62 */
/* bifurcation block 63 */
/* bifurcation block 64 */
/* bifurcation block 65 */
/* bifurcation block 66 */
/* bifurcation block 67 */
/* bifurcation block 68 */
/* bifurcation block 69 */
/* bifurcation block 70 */
/* bifurcation block 71 */
/* bifurcation block 72 */
/* bifurcation block 73 */
/* bifurcation block 74 */
/* bifurcation block 75 */
/* bifurcation block 76 */
/* bifurcation block 77 */
/* bifurcation block 78 */
/* bifurcation block 79 */
/* bifurcation block 80 */
/* bifurcation block 81 */
/* bifurcation block 82 */
/* bifurcation block 83 */
/* bifurcation block 84 */
/* bifurcation block 85 */
/* bifurcation block 86 */
/* bifurcation block 87 */
/* bifurcation block 88 */
/* bifurcation block 89 */
/* bifurcation block 90 */
/* bifurcation block 91 */
/* bifurcation block 92 */
/* bifurcation block 93 */
/* bifurcation block 94 */
/* bifurcation block 95 */
/* bifurcation block 96 */
/* bifurcation block 97 */
/* bifurcation block 98 */
/* bifurcation block 99 */
/* bifurcation block 100 */
/* bifurcation block 101 */
/* bifurcation block 102 */
/* bifurcation block 103 */
/* bifurcation block 104 */
/* bifurcation block 105 */
/* bifurcation block 106 */
/* bifurcation block 107 */
/* bifurcation block 108 */
/* bifurcation block 109 */
/* bifurcation block 110 */
/* bifurcation block 111 */
/* bifurcation block 112 */
/* bifurcation block 113 */
/* bifurcation block 114 */
/* bifurcation block 115 */
/* bifurcation block 116 */
/* bifurcation block 117 */
/* bifurcation block 118 */
/* bifurcation block 119 */
/* bifurcation block 120 */
/* bifurcation block 121 */
/* bifurcation block 122 */
/* bifurcation block 123 */
/* bifurcation block 124 */
/* bifurcation block 125 */
/* bifurcation block 126 */
/* bifurcation block 127 */
/* bifurcation block 128 */
/* bifurcation block 129 */
/* bifurcation block 130 */
/* bifurcation block 131 */
/* bifurcation block 132 */
/* bifurcation block 133 */
/* bifurcation block 134 */
/* bifurcation block 135 */
/* bifurcation block 136 */
/* bifurcation block 137 */
/* bifurcation block 138 */
/* bifurcation block 139 */
/* bifurcation block 140 */
/* bifurcation block 141 */
/* bifurcation block 142 */
/* bifurcation block 143 */
/* bifurcation block 144 */
/* bifurcation block 145 */
/* bifurcation block 146 */
/* bifurcation block 147 */
/* bifurcation block 148 */
/* bifurcation block 149 */
/* bifurcation block 150 */
/* bifurcation block 151 */
/* bifurcation block 152 */
/* bifurcation block 153 */
/* bifurcation block 154 */
/* bifurcation block 155 */
/* bifurcation block 156 */
/* bifurcation block 157 */
/* bifurcation block 158 */
/* bifurcation block 159 */
/* bifurcation block 160 */
/* bifurcation block 161 */
/* bifurcation block 162 */
/* bifurcation block 163 */
/* bifurcation block 164 */
/* bifurcation block 165 */
/* bifurcation block 166 */
/* bifurcation block 167 */
/* bifurcation block 168 */
/* bifurcation block 169 */
/* bifurcation block 170 */
/* bifurcation block 171 */
/* bifurcation block 172 */
/* bifurcation block 173 */
/* bifurcation block 174 */
/* bifurcation block 175 */
/* bifurcation block 176 */
/* bifurcation block 177 */
/* bifurcation block 178 */
/* bifurcation block 179 */
/* bifurcation block 180 */
/* bifurcation block 181 */
/* bifurcation block 182 */
/* bifurcation block 183 */
/* bifurcation block 184 */
/* bifurcation block 185 */
/* bifurcation block 186 */
/* bifurcation block 187 */
/* bifurcation block 188 */
/* bifurcation block 189 */
/* bifurcation block 190 */
/* bifurcation block 191 */
/* bifurcation block 192 */
/* bifurcation block 193 */
/* bifurcation block 194 */
/* bifurcation block 195 */
/* bifurcation block 196 */
/* bifurcation block 197 */
/* bifurcation block 198 */
/* bifurcation block 199 */
/* bifurcation block 200 */
double ncont_bifurcation_location_estimate(const NCONT_Bifurcation* b1, const NCONT_Bifurcation* b2) {
    if(!b1||!b2)return 0;return(b1->lambda+b2->lambda)/2.0;
}
int ncont_bifurcation_index(const NCONT_Bifurcation* b) {return b?b->point_index:-1;}
double ncont_bifurcation_get_lambda(const NCONT_Bifurcation* b) {return b?b->lambda:0;}
