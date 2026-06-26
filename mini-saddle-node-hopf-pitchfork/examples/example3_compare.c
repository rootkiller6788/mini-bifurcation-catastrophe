/* Compare all three codim-1 bifurcations via scanner */
#include "../include/bifurcation_detection.h"
#include <stdio.h>
#include <math.h>

static void compare_sys(const double* x, double* params, double* dx, int n) {
    /* System with tunable bifurcation: dx/dt = r*x + a - x^3 */
    (void)n;
    dx[0] = params[0] * x[0] + params[1] - x[0] * x[0] * x[0];
}

int main(void) {
    printf("=== Bifurcation Scanner: Compare All Types ===\n\n");

    /* Scan for saddle-node (vary a at r=0) */
    double x_sn[] = {1.0};
    double p_sn[] = {0.0, -0.5};  /* r=0, a varies */
    SaddleNodeBifurcation sn_res;
    int sn_rc = sn_detect(compare_sys, x_sn, p_sn, 1, 1, -2.0, 2.0, 0.1, 30, 1e-8, &sn_res);
    printf("Saddle-Node: %s at r=%.4f, x=%.4f\n",
           sn_rc >= 0 ? "FOUND" : "not found", sn_res.r_critical, sn_res.x_critical);

    /* Multi-type scanner */
    BifurcationScanner* bs = bscan_create();
    double x0[] = {0.5};
    double p0[] = {0.0, -0.3};
    int n_ev = bscan_detect_all(bs, compare_sys, x0, p0, 1, 0, -2.0, 2.0, 40, 30, 1e-6);

    printf("\nScanner found %d events across %d steps:\n", n_ev, bs->n_scan);
    for (int i = 0; i < bs->n; i++) {
        BifurcationEvent* ev = bscan_get_event(bs, i);
        printf("  Event %d: type=%s at param=%.4f\n", i, bscan_type_name(ev->type), ev->param);
    }

    int n_sn = bscan_count_type(bs, BIF_SN);
    int n_hf = bscan_count_type(bs, BIF_HOPF);
    int n_pf = bscan_count_type(bs, BIF_PITCHFORK);
    printf("\nCounts: SN=%d, Hopf=%d, Pitchfork=%d\n", n_sn, n_hf, n_pf);

    /* Check Z2 symmetry for pitchfork classification */
    bool z2 = pf_check_z2_symmetry(compare_sys, 1);
    printf("Z2 symmetry: %s\n", z2 ? "present" : "absent");

    bscan_print(bs);
    bscan_free(bs);
    printf("\nExample 3 PASSED\n");
    return 0;
}
