/* nft_app1.c -- L7: Nonlinear vibration normal form (Nayfeh 1973, Duffing).
 * Duffing oscillator: xddot + delta*xdot + alpha*x + beta*x^3 = F*cos(Omega*t)
 * Near primary resonance Omega ~ omega0, the method of multiple scales yields
 * a normal form for the slow amplitude dynamics which predicts:
 * - Backbone curve (nonlinear frequency-amplitude relation)
 * - Jump phenomenon (hysteresis in frequency response) */
#include "../include/nft_core.h"
#include <math.h>
#include <stdio.h>

int nft_duffing_backbone(double alpha, double beta, double delta, double F,
                         double* omega_vals, double* amp_vals, int n) {
    double omega0 = sqrt(alpha);
    printf("[nft_app1] Duffing oscillator: alpha=%.2f beta=%.2f delta=%.2f F=%.2f\n",
           alpha, beta, delta, F);
    printf("  Natural freq omega0=%.4f rad/s\n", omega0);

    /* Setup normal form matrix (linear part in action-angle) */
    NFTMatrix* A = nft_matrix_create(2, 2);
    A->data[0] = 0.0; A->data[1] = omega0;
    A->data[2] = -omega0; A->data[3] = 0.0;

    /* Normal form analysis */
    double ev_r[2], ev_i[2];
    int nev = nft_extract_eigenvalues(A, ev_r, ev_i);
    printf("  Eigenvalues: %.4f%+.4fi\n", ev_r[0], ev_i[0]);

    /* Compute backlash curve: Omega = omega0 + 3*beta/(8*omega0)*a^2 */
    for (int i = 0; i < n; i++) {
        double a = 0.1 + 2.9 * (double)i / (double)(n - 1);
        amp_vals[i] = a;
        /* Frequency response equation: [alpha - Omega^2 + 3/4*beta*a^2]^2 + (delta*Omega)^2 = (F/a)^2 */
        double denom = 1.5 * beta * a * a;
        omega_vals[i] = omega0 + denom / (2.0 * omega0);
    }

    /* Detect normal form bifurcation */
    NormalFormResult nf_res;
    nf_res.normal_form = A;
    Unfolding* uf = nft_unfolding_create(2);
    nft_compute_versal_unfolding(A, uf);
    BifurcationType bt = nft_classify_bifurcation(&nf_res, 0.0);
    printf("  Normal form bifurcation: %s\n", nft_bifurcation_name(bt));
    printf("  Tangent space dim: %d\n", nft_tangent_space_dimension(A, &nf_res));

    /* Output backbone curve */
    printf("  Backbone curve (frequency-amplitude):\n");
    for (int i = 0; i < n; i += n/6) {
        printf("    a=%.3f Omega=%.4f\n", amp_vals[i], omega_vals[i]);
    }

    nft_matrix_free(A);
    nft_unfolding_free(uf);
    return 0;
}

int main(void) {
    double om[30], amp[30];
    nft_duffing_backbone(1.0, 0.5, 0.05, 0.2, om, amp, 20);
    printf("=== nft_app1: Duffing Nonlinear Vibration ===\n");
    printf("Amplitude range: [%.2f, %.2f], freq range: [%.3f, %.3f]\n",
           amp[0], amp[19], om[0], om[19]);
    return 0;
}
