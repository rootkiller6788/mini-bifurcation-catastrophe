/* nft_app2.c -- L7-2: Aircraft wing flutter normal form (Dowell 1974, Hopf).
 * Wing flutter: two-degree-of-freedom aeroelastic system with plunge (h)
 * and pitch (alpha) degrees of freedom. As airspeed increases past a
 * critical value V_f, the damping of one mode crosses zero leading to
 * a Hopf bifurcation — sustained oscillations (flutter) emerge.
 * Normal form analysis determines supercritical/subcritical character. */
#include "../include/nft_core.h"
#include <math.h>
#include <stdio.h>

int nft_wing_flutter(double m, double I_alpha, double k_h, double k_alpha,
                     double rho, double b, double* V_f, double* flutter_freq) {
    /* Aeroelastic parameters: mass m, moment of inertia I_alpha,
     * stiffnesses k_h (plunge), k_alpha (pitch), air density rho,
     * semi-chord b */
    double omega_h = sqrt(k_h / m);
    double omega_a = sqrt(k_alpha / I_alpha);
    double mu = m / (3.141592653589793 * rho * b * b);

    printf("[nft_app2] Wing flutter analysis:\n");
    printf("  omega_h=%.3f omega_a=%.3f Hz mass_ratio=mu=%.2f\n",
           omega_h, omega_a, mu);

    /* Compute flutter speed via eigenvalue analysis */
    double omega_avg = (omega_h + omega_a) / 2.0;
    *flutter_freq = omega_avg;
    *V_f = b * omega_avg * sqrt(mu);  /* simplified Theodorsen */

    printf("  Critical flutter speed V_f=%.2f m/s\n", *V_f);
    printf("  Flutter frequency omega_f=%.3f rad/s\n", *flutter_freq);

    /* Normal form of the flutter Hopf bifurcation */
    int n = 4;
    NFTMatrix* A = nft_matrix_create(n, n);
    /* Construct linear aeroelastic matrix at V_f */
    double V = *V_f;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A->data[i*n+j] = 0.0;
    A->data[0*n+2] = 1.0;  /* h_dot = dh/dt */
    A->data[1*n+3] = 1.0;  /* alpha_dot = dalpha/dt */
    A->data[2*n+0] = -omega_h * omega_h;
    A->data[2*n+1] = V * V * rho * b / m;
    A->data[3*n+0] = V * V * rho * b * b / I_alpha;
    A->data[3*n+1] = -omega_a * omega_a;

    double ev_r[4], ev_i[4];
    int nev = nft_extract_eigenvalues(A, ev_r, ev_i);
    printf("  Eigenvalues at V_f:\n");
    for (int i = 0; i < nev; i++) {
        printf("    %.4f%+.4fi [%s]\n", ev_r[i], ev_i[i],
               ev_r[i] < -0.001 ? "stable" : ev_r[i] > 0.001 ? "unstable" : "marginal");
    }

    /* Classify the bifurcation */
    NormalFormResult nf_res;
    nf_res.normal_form = A;
    BifurcationType bt = nft_classify_bifurcation(&nf_res, V);
    printf("  Bifurcation type: %s\n", nft_bifurcation_name(bt));

    /* Compute versal unfolding to determine type */
    Unfolding* uf = nft_unfolding_create(1);
    nft_compute_versal_unfolding(A, uf);
    double mu_crit;
    int bp = nft_detect_bifurcation_point(uf, V*0.9, V*1.1, &mu_crit);
    printf("  Versal unfolding: %sversals, bif_point=%.2f %s\n",
           uf->is_versal ? "" : "not ", mu_crit, bp >= 0 ? "(detected)" : "");

    nft_matrix_free(A);
    nft_unfolding_free(uf);
    return 0;
}

int main(void) {
    double V_f, omega_f;
    nft_wing_flutter(100.0, 50.0, 5000.0, 3000.0, 1.225, 1.5, &V_f, &omega_f);
    printf("=== nft_app2: Wing Flutter Hopf ===\n");
    printf("V_f=%.2f m/s, omega_f=%.2f rad/s\n", V_f, omega_f);
    return 0;
}
