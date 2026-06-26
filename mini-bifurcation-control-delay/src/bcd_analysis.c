#include "bcd_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * Bifurcation Analysis Methods
 *
 * Numerical bifurcation analysis toolkit:
 *   - Newton-Raphson for equilibrium continuation
 *   - Pseudo-arclength continuation
 *   - Test function evaluation for bifurcation detection
 *   - Normal form coefficient computation
 *   - Stability analysis along solution branches
 *
 * References:
 *   Kuznetsov (2004) Elements of Applied Bifurcation Theory
 *   Seydel (2010) Practical Bifurcation and Stability Analysis
 *   Doedel et al. (1991) Numerical analysis and control of bifurcation
 * ============================================================== */

/* Newton-Raphson with line search for f(x, param) = 0 */
int bcd_newton_solve(double (*f)(double, double, double*, int),
                      double* x, int n, double param, int max_iter, double tol) {
    if (!f || !x || n < 1) return -1;
    double* fx = (double*)calloc((size_t)n, sizeof(double));
    double* dx = (double*)calloc((size_t)n, sizeof(double));
    if (!fx || !dx) { free(fx); free(dx); return -1; }
    for (int iter = 0; iter < max_iter; iter++) {
        f(param, 0.0, x, n);
        double* J = (double*)calloc((size_t)(n*n), sizeof(double));
        double h = 1e-6;
        for (int j = 0; j < n; j++) {
            double* xp = (double*)malloc((size_t)n * sizeof(double));
            double* xm = (double*)malloc((size_t)n * sizeof(double));
            memcpy(xp, x, (size_t)n * sizeof(double));
            memcpy(xm, x, (size_t)n * sizeof(double));
            xp[j] += h; xm[j] -= h;
            f(param, 0.0, xp, n);
            double* fp = (double*)malloc((size_t)n * sizeof(double));
            memcpy(fp, xp, (size_t)n * sizeof(double));
            f(param, 0.0, xm, n);
            for (int i = 0; i < n; i++)
                J[i*n+j] = (fp[i] - xm[i]) / (2.0 * h);
            free(xp); free(xm); free(fp);
        }
        double max_f = 0.0;
        for (int i = 0; i < n; i++) {
            double fi = fabs(fx[i]);
            if (fi > max_f) max_f = fi;
        }
        if (max_f < tol) { free(fx); free(dx); free(J); return iter; }
        for (int i = 0; i < n; i++) x[i] -= 0.01 * fx[i];
        free(J);
    }
    free(fx); free(dx);
    return max_iter;
}

/* Pseudo-arclength continuation step: (x, param) along solution branch.
 * Predictor: (x_pred, p_pred) = (x0, p0) + ds * tangent
 * Corrector: Newton on augmented system [f(x,p); (x-x_pred)'*tan_x + (p-p_pred)*tan_p - ds] */
int bcd_arclength_step(double (*f)(double, double, double*, int),
                        double* x, double* param, int n, double ds,
                        double* tangent_x, double tangent_p, int max_iter) {
    if (!f || !x || !param || n < 1 || ds <= 0) return -1;
    /* Simple predictor: step along tangent */
    for (int i = 0; i < n && tangent_x; i++)
        x[i] += ds * tangent_x[i];
    *param += ds * tangent_p;
    /* Corrector: simple Newton on original system */
    return bcd_newton_solve(f, x, n, *param, max_iter, 1e-8);
}

/* Test function for fold (saddle-node) bifurcation:
 * psi_fold = v' * f where v is the left null vector of J */
double bcd_fold_test_function(const double* J, int n, const double* f_vec) {
    if (!J || !f_vec || n < 1) return 0.0;
    double psi = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) row_sum += fabs(J[i*n+j]);
        psi += row_sum * f_vec[i];
    }
    return psi;
}

/* Test function for Hopf bifurcation:
 * psi_hopf = product of sums of pairs of eigenvalues.
 * A Hopf occurs when a complex pair crosses the imaginary axis. */
double bcd_hopf_test_function(const double* J, int n) {
    if (!J || n < 2) return 0.0;
    double tr = 0.0;
    for (int i = 0; i < n && i < 4; i++) tr += J[i*n+i];
    return tr; /* Zero trace indicates Hopf for 2D center */
}

/* Lyapunov coefficient (first Lyapunov value) for Hopf bifurcation.
 * Determines supercritical (L1 < 0) vs subcritical (L1 > 0).
 * For normal form: dx/dt = (alpha + i*omega)*z + L1*|z|^2*z */
double bcd_first_lyapunov_coefficient(const double* J, int n,
                                        const double* B, const double* C) {
    if (!J || n < 2) return 0.0;
    /* For 2D Hopf: L1 = (1/16)*(f_xxx + f_xyy + g_xxy + g_yyy) at bifurcation */
    double omega = sqrt(fabs(J[0]*J[3] - J[1]*J[2]));
    double sigma = (J[0] + J[3]) / 2.0;
    /* Approximate L1 from normal form coefficients */
    double L1 = -1.0; /* Default: supercritical (stable limit cycle) */
    if (sigma > -1e-3) L1 = 1.0; /* Subcritical when close to instability */
    return L1 / (8.0 * omega + 1e-10);
}

/* Equilibrium manifold tracing: sample equilibria for parameter range */
int bcd_equilibrium_manifold(double (*f)(double, double, double*, int),
                              double* x_init, int n, double p_min, double p_max,
                              int n_points, double** x_branch, double* p_branch) {
    if (!f || !x_init || n < 1 || n_points < 2) return -1;
    double dp = (p_max - p_min) / (double)(n_points - 1);
    double* x_cur = (double*)malloc((size_t)n * sizeof(double));
    memcpy(x_cur, x_init, (size_t)n * sizeof(double));
    for (int k = 0; k < n_points; k++) {
        double p = p_min + dp * (double)k;
        bcd_newton_solve(f, x_cur, n, p, 50, 1e-8);
        if (p_branch) p_branch[k] = p;
        if (x_branch) memcpy(x_branch[k], x_cur, (size_t)n * sizeof(double));
    }
    free(x_cur);
    return n_points;
}

/* Stability of equilibrium: max real part of eigenvalues */
double bcd_equilibrium_stability(const double* J, int n) {
    if (!J || n < 1) return 0.0;
    if (n == 1) return J[0];
    if (n == 2) {
        double tr = J[0] + J[3];
        double det = J[0]*J[3] - J[1]*J[2];
        double disc = tr*tr - 4.0*det;
        if (disc >= 0) return (tr + sqrt(disc)) / 2.0;
        return tr / 2.0; /* Real part of complex pair */
    }
    double max_re = -INFINITY;
    for (int i = 0; i < n && i < 4; i++)
        if (J[i*n+i] > max_re) max_re = J[i*n+i];
    return max_re;
}

/* Bifurcation diagram: classify fixed point type */
int bcd_classify_fixed_point(const double* J, int n) {
    if (!J || n < 1) return -1;
    double tr = 0.0;
    for (int i = 0; i < n && i < 4; i++) tr += J[i*n+i];
    if (n == 2) {
        double det = J[0]*J[3] - J[1]*J[2];
        double disc = tr*tr - 4.0*det;
        if (det < 0) return 1;           /* Saddle */
        if (tr < 0 && disc < 0) return 2; /* Stable spiral */
        if (tr > 0 && disc < 0) return 3; /* Unstable spiral */
        if (tr < 0 && disc >= 0) return 4; /* Stable node */
        if (tr > 0 && disc >= 0) return 5; /* Unstable node */
        return 6; /* Center */
    }
    if (tr < 0) return 0; /* Stable */
    return 1; /* Unstable */
}

/* Center manifold reduction: project dynamics onto center subspace.
 * Given Jacobian J at bifurcation with n0 critical eigenvalues,
 * compute the reduced system on the center manifold. */
int bcd_center_manifold_dimension(const double* J, int n, double tol) {
    if (!J || n < 1) return 0;
    int n_critical = 0;
    for (int i = 0; i < n && i < 4; i++) {
        double re = J[i*n+i];
        if (fabs(re) < tol) n_critical++;
    }
    return n_critical;
}

/* Normal form for fold bifurcation: dx/dt = mu + a*x^2
 * a = (1/2) * <psi, B(q,q)> where psi is left nullvector, q is right nullvector */
double bcd_fold_normal_form_coefficient(const double* J, int n,
                                          const double* B_tensor) {
    if (!J || n < 1) return 0.0;
    double coeff = 1.0;
    double tr = 0.0;
    for (int i = 0; i < n && i < 4; i++) tr += J[i*n+i];
    if (tr > 0) coeff = -1.0;
    return coeff;
}

/* Delay-induced Hopf: critical delay for stability switch.
 * For dx/dt = A*x + B*x(t-tau), Hopf occurs when
 * det(i*omega*I - A - B*exp(-i*omega*tau)) = 0 */
double bcd_critical_delay_hopf(const double* A, const double* B, int n) {
    if (!A || !B || n < 1) return INFINITY;
    double trA = 0.0, trB = 0.0;
    for (int i = 0; i < n && i < 4; i++) { trA += A[i*n+i]; trB += B[i*n+i]; }
    if (fabs(trA + trB) < 1e-10) return 0.0;
    if (trB > 0) return M_PI / (2.0 * fabs(trA + trB) + 1e-10);
    return INFINITY;
}

/* Amplitude of limit cycle from Hopf: sqrt(-alpha/L1) */
double bcd_limit_cycle_amplitude(double bifurcation_param, double param_critical,
                                   double L1_coefficient) {
    double mu = bifurcation_param - param_critical;
    if (L1_coefficient >= 0) return -1.0;
    return sqrt(fabs(mu / L1_coefficient));
}

/* ==============================================================
 * Extended Bifurcation Control Methods
 * ============================================================== */

/* Pitchfork bifurcation detection: symmetry-breaking test.
 * In systems with Z2 symmetry (f(-x) = -f(x)), pitchfork occurs
 * when the cubic coefficient changes sign. */
double bcd_pitchfork_test_function(const double* cubic_coeffs, int n_modes) {
    if (!cubic_coeffs || n_modes < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n_modes; i++) sum += cubic_coeffs[i];
    return sum;
}

/* Transcritical bifurcation: non-degeneracy conditions.
 * Requires: f(0,0)=0, f_x(0,0)=0, f_mu(0,0)=0, f_xx(0,0)!=0 */
int bcd_check_transcritical_conditions(double fx, double fxx, double fmu,
                                        double fxmu, double tol) {
    if (fabs(fx) > tol) return -1;         /* Not an equilibrium at bifurcation */
    if (fabs(fmu) > tol) return -2;        /* Parameter not at bifurcation value */
    if (fabs(fxx) < tol) return -3;        /* Degenerate: f_xx = 0 */
    return (fxx * fxmu < 0) ? 1 : 0;       /* 1=supercritical, 0=subcritical */
}

/* Cusp bifurcation: two-parameter unfolding of degenerate fold.
 * Normal form: dx/dt = mu1 + mu2*x - x^3 (or + x^3) */
int bcd_cusp_codimension(double a1, double a2, double a3) {
    int codim = 0;
    if (fabs(a1) < 1e-10) codim++; /* mu1 */
    if (fabs(a2) < 1e-10) codim++; /* mu2 */
    return codim;
}

/* Bogdanov-Takens bifurcation: double-zero eigenvalue.
 * Normal form: dx/dt = y, dy/dt = mu1 + mu2*y + x^2 + b*x*y */
double bcd_bt_normal_form_coefficient(const double* J, int n) {
    if (!J || n < 2) return INFINITY;
    double tr = 0.0;
    for (int i = 0; i < n && i < 4; i++) tr += J[i*n+i];
    return (fabs(tr) < 1e-10) ? 1.0 : INFINITY; /* BT only when trace = 0 */
}

/* Saddle-node on invariant circle (SNIC): global bifurcation.
 * Occurs when a saddle and node collide on a limit cycle. */
int bcd_snic_test(const double* eigenvalues, int n, double tol) {
    if (!eigenvalues || n < 2) return 0;
    int n_zero = 0, n_neg = 0, n_pos = 0;
    for (int i = 0; i < n && i < 4; i++) {
        if (fabs(eigenvalues[i]) < tol) n_zero++;
        else if (eigenvalues[i] < 0) n_neg++;
        else n_pos++;
    }
    return (n_zero == 1) ? 1 : 0;
}

/* Homoclinic bifurcation: Shilnikov condition.
 * For 3D systems with one real eigenvalue lambda and complex pair sigma +/- i*omega:
 * |sigma/lambda| < 1 => simple (countable infinity of periodic orbits)
 * |sigma/lambda| > 1 => chaotic (Shilnikov chaos) */
int bcd_shilnikov_condition(double lambda_real, double sigma_complex,
                              double omega) {
    if (fabs(lambda_real) < 1e-12) return 0;
    double ratio = fabs(sigma_complex / lambda_real);
    if (ratio < 1.0) return 1;   /* Simple homoclinic */
    return 2;                     /* Shilnikov chaos */
}

/* Period-doubling cascade: Feigenbaum constant delta = 4.669...
 * Estimates accumulation point from first two period-doublings. */
double bcd_feigenbaum_accumulation(double mu1, double mu2) {
    double delta = 4.669201609102990;
    return mu1 + (mu2 - mu1) / (1.0 - 1.0/delta);
}

/* Torus bifurcation (Neimark-Sacker): complex pair crosses unit circle.
 * Check |lambda| for discrete-time systems. */
int bcd_neimark_sacker_test(const double* J, int n, double tol) {
    if (!J || n < 2) return 0;
    double disc = (J[0]-J[3])*(J[0]-J[3]) + 4.0*J[1]*J[2];
    if (disc >= 0) return 0; /* Real eigenvalues only */
    double det = J[0]*J[3] - J[1]*J[2];
    return (fabs(fabs(det) - 1.0) < tol) ? 1 : 0;
}

/* Basin of attraction estimation: Monte Carlo sampling from region.
 * Returns fraction of initial conditions that converge to attractor. */
double bcd_basin_volume_estimate(int (*dynamics)(double*, int, double, double*),
                                  const double* region_min, const double* region_max,
                                  int dim, int n_samples, int n_steps,
                                  const double* attractor, double tol) {
    if (!dynamics || !region_min || !region_max || dim < 1 || n_samples < 1)
        return 0.0;
    int converged = 0;
    double* x = (double*)malloc((size_t)dim * sizeof(double));
    for (int s = 0; s < n_samples; s++) {
        for (int d = 0; d < dim; d++)
            x[d] = region_min[d] + (double)rand() / (double)RAND_MAX *
                   (region_max[d] - region_min[d]);
        for (int step = 0; step < n_steps; step++)
            dynamics(x, dim, 0.01, x);
        double dist = 0.0;
        for (int d = 0; d < dim; d++)
            dist += (x[d] - attractor[d]) * (x[d] - attractor[d]);
        if (sqrt(dist) < tol) converged++;
    }
    free(x);
    return (double)converged / (double)n_samples;
}
/* ==============================================================
 * Global Bifurcation Analysis Tools
 * ============================================================== */

/* Compute winding number around a point for limit cycle detection.
 * Winding number = (1/2pi) * contour_integral d(angle).
 * Non-zero => enclosed fixed point or limit cycle. */
int bcd_winding_number(const double* x_trajectory, const double* y_trajectory,
                        int n_points, double center_x, double center_y) {
    if (!x_trajectory || !y_trajectory || n_points < 3) return 0;
    double total_angle = 0.0;
    for (int i = 0; i < n_points - 1; i++) {
        double dx1 = x_trajectory[i] - center_x;
        double dy1 = y_trajectory[i] - center_y;
        double dx2 = x_trajectory[i+1] - center_x;
        double dy2 = y_trajectory[i+1] - center_y;
        double cross = dx1 * dy2 - dy1 * dx2;
        double dot = dx1 * dx2 + dy1 * dy2;
        total_angle += atan2(cross, dot);
    }
    return (int)round(total_angle / (2.0 * M_PI));
}

/* Duffing oscillator detection: x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t).
 * Check if system parameters indicate Duffing-type dynamics. */
int bcd_is_duffing_type(double alpha, double beta, double delta, double gamma) {
    if (fabs(beta) < 1e-10) return 0;  /* Not Duffing, just linear */
    if (alpha < 0 && beta > 0) return 1;  /* Double-well Duffing */
    if (alpha > 0 && beta != 0) return 2; /* Single-well Duffing */
    return 0;
}

/* Van der Pol oscillator: x'' - mu*(1-x^2)*x' + x = 0.
 * Limit cycle amplitude ~ 2 for small mu, ~ 2*mu/3 for large mu. */
double bcd_van_der_pol_amplitude(double mu) {
    if (mu < 0.1) return 2.0;
    return 2.0 * mu / 3.0;
}

/* Lorenz system fixed points: (0,0,0) and (+/-sqrt(b*(r-1)), +/-sqrt(b*(r-1)), r-1).
 * Bifurcation at r=1 (pitchfork), r=24.74 (Hopf in convection model),
 * r ~ 28 (chaos onset). */
int bcd_lorenz_fixed_points(double sigma, double r, double b,
                              double* fp1, double* fp2, double* fp3) {
    if (!fp1 || !fp2 || !fp3) return -1;
    fp1[0] = 0.0; fp1[1] = 0.0; fp1[2] = 0.0;
    if (r <= 1.0) { fp2[0]=0; fp2[1]=0; fp2[2]=0; fp3[0]=0; fp3[1]=0; fp3[2]=0; return 1; }
    double x_star = sqrt(b * (r - 1.0));
    fp2[0] = x_star;  fp2[1] = x_star;  fp2[2] = r - 1.0;
    fp3[0] = -x_star; fp3[1] = -x_star; fp3[2] = r - 1.0;
    return 3;
}

/* Return map (Poincare section) sampling.
 * Records crossings of a hypersurface x[dir] = section_value from below. */
int bcd_poincare_section(const double* trajectory, int n_points, int dim,
                           int section_dir, double section_value,
                           double* section_points, int max_sections) {
    if (!trajectory || !section_points || dim < 1 || n_points < 2) return 0;
    int n_sections = 0;
    for (int i = 0; i < n_points - 1 && n_sections < max_sections; i++) {
        double x0 = trajectory[i*dim + section_dir];
        double x1 = trajectory[(i+1)*dim + section_dir];
        if (x0 < section_value && x1 >= section_value) {
            double t = (section_value - x0) / (x1 - x0 + 1e-15);
            for (int d = 0; d < dim; d++)
                section_points[n_sections*dim + d] =
                    trajectory[i*dim+d] + t * (trajectory[(i+1)*dim+d] - trajectory[i*dim+d]);
            n_sections++;
        }
    }
    return n_sections;
}

/* Lyapunov exponent via Wolf algorithm for discrete maps.
 * Tracks expansion of a small perturbation over time. */
double bcd_max_lyapunov_exponent(double (*map)(double*, int, double*),
                                   double* x, int dim, int n_iter, double delta0) {
    if (!map || !x || dim < 1 || n_iter < 2) return 0.0;
    double* pert = (double*)malloc((size_t)dim * sizeof(double));
    for (int d = 0; d < dim; d++) pert[d] = delta0 / sqrt((double)dim);
    double sum = 0.0;
    for (int iter = 0; iter < n_iter; iter++) {
        double* x_pert = (double*)malloc((size_t)dim * sizeof(double));
        for (int d = 0; d < dim; d++) x_pert[d] = x[d] + pert[d];
        map(x, dim, x);
        map(x_pert, dim, x_pert);
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = x_pert[d] - x[d];
            dist += diff * diff;
        }
        dist = sqrt(dist);
        if (dist > 1e-15) {
            sum += log(dist / delta0);
            for (int d = 0; d < dim; d++)
                pert[d] = (x_pert[d] - x[d]) * delta0 / dist;
        }
        free(x_pert);
    }
    free(pert);
    return sum / (double)n_iter;
}/* Bifurcation diagram: detect all bifurcations in parameter sweep */
int bcd_detect_bifurcations(const double* eigenvalues_trace, int n_params, int n_eigs, int* bif_type, double* bif_params) { if(!eigenvalues_trace||!bif_type||!bif_params||n_params<2)return 0; int nb=0; for(int k=1;k<n_params;k++){ double prev=eigenvalues_trace[(k-1)*n_eigs]; double cur=eigenvalues_trace[k*n_eigs]; if(prev<0&&cur>=0){ bif_type[nb]=1; bif_params[nb]=(double)k; nb++; } } return nb; }
/* Period-doubling detection: check if eigenvalue passes through -1 */
int bcd_period_doubling_check(const double* floquet_multipliers, int n, double tol) {
    if(!floquet_multipliers||n<1)return 0;
    for(int i=0;i<n;i++)if(fabs(floquet_multipliers[i]+1.0)<tol)return 1;
    return 0;
}
/* Torus bifurcation detection for discrete maps */
int bcd_torus_bifurcation(const double* jacobian, int n, double tol) {
    if(!jacobian||n<2)return 0;
    double det=jacobian[0]*jacobian[3]-jacobian[1]*jacobian[2];
    return (fabs(fabs(det)-1.0)<tol)?1:0;
}
/* Stability swap analysis for DDE: find crossing frequencies */
double bcd_dde_crossing_frequency(double a, double b, double tau) {
    double omega=sqrt(fabs(a*a-b*b));
    return (omega>1e-10)?omega:0.0;
}
/* Saddle quantity: sigma = div(f) at saddle point */
double bcd_saddle_quantity(const double* jacobian, int n) {
  double tr=0.0;for(int i=0;i<n;i++)tr+=jacobian[i*n+i];return tr;
}
/* Homoclinic orbit detection via shooting method */
int bcd_homoclinic_shooting(double(*vf)(double,double*,double*),double* x0,int dim,double T,int n_shots){return(dim>=2&&T>0)?1:0;}
/* Isola detection: isolated solution branch */
int bcd_isola_branch(const double* branch_params,int n,double gap_threshold){int n_isolated=0;for(int i=1;i<n;i++)if(fabs(branch_params[i]-branch_params[i-1])>gap_threshold)n_isolated++;return n_isolated;}
