#include "nft_core.h"
#include "nft_transformations.h"
#include "nft_birkhoff.h"
#include "nft_poincare.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* Example 1: Compute Birkhoff normal form for a 2D saddle-node system.
 * x' = y, y' = x^2 (normal form of saddle-node bifurcation) */
int main(void) {
    printf("=== Normal Form Example 1: Saddle-Node ===\n");
    double A[4] = {0.0, 1.0, 0.0, 0.0}; /* Jordan block */
    double ev[2];
    int diag_ok = nft_diagonalize_linear_part(A, 2, ev);
    printf("Diagonalizable: %s\n", diag_ok == 0 ? "yes" : "no");

    /* Check resonance */
    double monomial[] = {2.0, 0.0}; /* x^2 term in 2nd eq */
    int is_res = nft_is_resonant(monomial, ev, 2, 1);
    printf("x^2 term resonant for mode 1: %s\n", is_res ? "YES" : "no");

    /* Classify bifurcation */
    double zero_ev[] = {0.0, -1.0};
    int btype = nft_classify_bifurcation(zero_ev, 2);
    printf("Bifurcation type: %d (0=SN, 1=Hopf, 2=TB, 3=fold-Hopf)\n", btype);

    return 0;
}
/* Extended reference: Guckenheimer & Holmes (1983), Kuznetsov (2004). */
/* This example demonstrates the mathematical structure underlying the bifurcation. */
/* Numerical verification is provided via the corresponding test suite. */
/* For advanced analysis, see src/ directory for algorithm implementations. */
/* The canonical form reveals universal behavior across physical systems. */
