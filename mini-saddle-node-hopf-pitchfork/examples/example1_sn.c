/* Saddle-node bifurcation: detect and analyze dx/dt = r + x^2 */
#include "../include/saddle_node.h"
#include <stdio.h>
#include <math.h>

static void sn_system(const double* x, double* params, double* dx, int n) {
    (void)n;
    dx[0] = params[0] + x[0] * x[0];  /* dx/dt = r + x^2 */
}

int main(void) {
    printf("=== Saddle-Node Bifurcation Detection ===\n\n");

    double x0[] = {0.5};
    double params[] = {-1.0, 0.0};
    SaddleNodeBifurcation* sn = sn_create();

    int rc = sn_detect(sn_system, x0, params, 1, 0, -1.5, 1.5, 0.1, 30, 1e-8, sn);
    if (rc < 0) { printf("Detection FAILED\n"); sn_free(sn); return 1; }

    printf("Detected at r_critical = %.6f, x_critical = %.6f\n", sn->r_critical, sn->x_critical);
    printf("Type: %s\n", sn->is_fold ? "Fold (saddle-node)" : "Other");

    /* Verify: compute branches around bifurcation */
    double xb[] = {0.5};
    double pb[] = {-0.5, 0.0};
    int npts = sn_compute_branches(sn, sn_system, xb, pb, 1, 0, -2.0, 2.0, 40, 30, 1e-6);
    printf("Branches computed: %d points\n", npts);

    if (npts > 0 && sn->stable_branch && sn->unstable_branch) {
        printf("\nBranch data (first 5 points):\n");
        for (int i = 0; i < npts && i < 5; i++) {
            printf("  p[%d]=%.3f  stable=%.4f  unstable=%.4f\n",
                   i, sn->param_range[i], sn->stable_branch[i], sn->unstable_branch[i]);
        }
    }

    /* Check eigenvalue condition */
    Matrix* J = mat_create(1, 1);
    mat_set(J, 0, 0, 2.0 * sn->x_critical);
    double fold_cond = sn_fold_condition(J);
    printf("\nFold condition det(J) at crit pt: %.6f (expect near 0)\n", fold_cond);
    mat_free(J);

    sn_print(sn);
    sn_free(sn);
    printf("\nExample 1 PASSED\n");
    return 0;
}
