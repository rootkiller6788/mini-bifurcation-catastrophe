#ifndef NCONT_BRANCH_H
#define NCONT_BRANCH_H
#include "ncont_core.h"
#include "ncont_bifurcation.h"

/* ==============================================================
 * ncont_branch.h ? Branch Switching at Bifurcation Points
 *
 * At a bifurcation point, multiple solution branches intersect.
 * Branch switching finds initial guesses for the new branch.
 *
 * Methods:
 *   - Eigenvector perturbation (algebraic bifurcation equation)
 *   - Lyapunov-Schmidt reduction
 *   - Perturbed bifurcation parameter
 *
 * References:
 *   Keller (1977)
 *   Seydel (2010) Ch.5
 *   Kuznetsov (2004) Ch.10
 * ============================================================== */

/* Branch switching result */
typedef struct {
    double x_new[NCONT_MAX_DIM];
    int n;
    double lambda_new;
    NCONT_BifurcationType bif_type;
    bool success;
    int direction;            /* +1 or -1 for which branch direction */
} NCONT_BranchSwitch;

/* Switch branches at a limit point (fold) */
NCONT_BranchSwitch ncont_switch_limit_point(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params, double ds);

/* Switch branches at a branch point (pitchfork/transcritical) */
NCONT_BranchSwitch ncont_switch_branch_point(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params, double ds);

/* Switch branches at a Hopf point (to periodic orbit branch) */
NCONT_BranchSwitch ncont_switch_hopf(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Point* bif_point, int n, void* params,
    double omega_guess, double ds);

/* Generic branch switching based on bifurcation type */
NCONT_BranchSwitch ncont_switch_branch(NCONT_System F, NCONT_Jacobian J,
    const NCONT_Bifurcation* bif, const NCONT_Point* bif_point,
    int n, void* params, double ds);

/* Compute null vector of Jacobian at bifurcation */
bool ncont_compute_null_vector(NCONT_Jacobian J, const double* x, int n,
    double lambda, void* params, double* null_vector);

/* Perturb solution along eigenvector direction */
void ncont_perturb_along_vector(const double* x, int n,
    const double* direction, double epsilon, double* x_new);

/* Check if branch switch was successful */
bool ncont_branch_switch_success(const NCONT_BranchSwitch* bs);

#endif

























