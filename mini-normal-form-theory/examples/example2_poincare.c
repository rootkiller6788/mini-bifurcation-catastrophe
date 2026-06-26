#include "nft_core.h"
#include "nft_transformations.h"
#include "nft_poincare.h"
#include <stdio.h>
#include <math.h>

/* Example 2: Poincare-Dulac theorem demo.
 * Reduce a cubic nonlinearity to its normal form by
 * eliminating non-resonant terms via near-identity transform. */
int main(void) {
    printf("=== Normal Form Example 2: Poincare-Dulac ===\n");

    /* Linear part: A = diag(1, -2) */
    double A[4] = {1.0, 0.0, 0.0, -2.0};
    double ev[2];
    nft_diagonalize_linear_part(A, 2, ev);
    printf("Eigenvalues: %.2f, %.2f\n", ev[0], ev[1]);

    /* Check which cubic monomials are resonant */
    double mono1[] = {2.0, 1.0}; /* x^2*y in eq 1 */
    printf("x^2*y in eq1 resonant: %s\n",
           nft_is_resonant(mono1, ev, 2, 0) ? "YES" : "no");

    double mono2[] = {1.0, 0.0}; /* x in eq 2 */
    printf("x in eq2 resonant: %s\n",
           nft_is_resonant(mono2, ev, 2, 1) ? "YES" : "no");

    /* Normal form coefficients */
    double coeff_arr[6] = {0.1, -0.3, 0.5, 0.2, -0.1, 0.4};
    double* coeffs[4];
    for (int i = 0; i < 4; i++) coeffs[i] = &coeff_arr[i];
    double nf_out[4] = {0};
    double* nf[4];
    for (int i = 0; i < 4; i++) nf[i] = &nf_out[i];
    nft_compute_normal_form(A, 2, (const double**)coeffs, 3, nf);
    printf("Normal form coeffs: [%.3f, %.3f, %.3f, %.3f]\n",
           nf_out[0], nf_out[1], nf_out[2], nf_out[3]);

    return 0;
}
