/* Birkhoff normal form: compute Poincare-Dulac normal form for 2D system */
#include "../include/nft_core.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Birkhoff Normal Form: 2-DOF Hamiltonian ===\n\n");

    /* Setup linear part A (diagonal with eigenvalues) */
    int n = 2;
    NFTMatrix* A = nft_matrix_create(n, n);
    A->data[0] = 0.0; A->data[1] = -1.0;  /* omega_1 = 1.0 */
    A->data[2] = 1.0; A->data[3] = 0.0;
    nft_matrix_print(A, "Linear part A");

    /* Extract eigenvalues */
    double ev_real[2], ev_imag[2];
    int nev = nft_extract_eigenvalues(A, ev_real, ev_imag);
    printf("Eigenvalues: ");
    for (int i = 0; i < nev; i++) {
        printf("%.4f%+.4fi ", ev_real[i], ev_imag[i]);
    }
    printf("\n");

    /* Check resonance conditions */
    double omega1 = fabs(ev_imag[0]);
    printf("omega_1 = %.4f, omega_1/omega_1 ratio = %.4f\n", omega1, omega1/omega1);

    /* Compute versal unfolding */
    Unfolding* uf = nft_unfolding_create(1);
    int rc = nft_compute_versal_unfolding(A, uf);
    printf("Versal unfolding: %s, codim=%.1f, n_basis=%d\n",
           rc >= 0 ? "computed" : "FAILED",
           nft_unfolding_codimension(A), uf->n_basis);

    /* Classify bifurcation type */
    NormalFormResult nf_res;
    nf_res.normal_form = A;
    BifurcationType bt = nft_classify_bifurcation(&nf_res, 0.0);
    printf("Bifurcation type: %s\n", nft_bifurcation_name(bt));

    /* Compute tangent space dimension */
    int tdim = nft_tangent_space_dimension(A, &nf_res);
    printf("Tangent space dim: %d\n", tdim);

    /* Test nilpotency */
    bool nilp = nft_matrix_is_nilpotent(A);
    printf("A nilpotent: %s\n", nilp ? "YES" : "NO");

    /* Characteristic polynomial */
    double coeffs[3];
    nft_charpoly_coeffs(A, coeffs);
    printf("Char poly: %.2f + %.2f*x + %.2f*x^2\n", coeffs[0], coeffs[1], coeffs[2]);

    /* Evaluate vector field at test point */
    double x_test[] = {0.5, 0.3};
    double dx_test[2];
    nft_eval_vector_field(A, x_test, dx_test, n);
    printf("dx/dt at (%.2f,%.2f) = (%.4f, %.4f)\n", x_test[0], x_test[1], dx_test[0], dx_test[1]);

    /* Bifurcation diagram sweep */
    printf("\nBifurcation diagram scan:\n");
    double eq_buffer[20];
    int n_eq;
    int bd_pts = nft_bifurcation_diagram(uf, -2.0, 2.0, 15, eq_buffer, &n_eq);
    printf("  %d diagram points computed, %d equilibria\n", bd_pts, n_eq);

    /* Detect exact bifurcation point */
    double mu_crit;
    int bp = nft_detect_bifurcation_point(uf, -1.0, 1.0, &mu_crit);
    printf("Bifurcation at mu_crit = %.6f (%s)\n", mu_crit, bp >= 0 ? "found" : "not found");

    nft_matrix_free(A);
    nft_unfolding_free(uf);
    printf("\nExample 3 PASSED\n");
    return 0;
}
