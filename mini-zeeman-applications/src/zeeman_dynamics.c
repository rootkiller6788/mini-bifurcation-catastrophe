#include "zeeman_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Jump Detection
 * ============================================================== */

double zeeman_find_delay_jump(double a_start, double a_end, double b_fixed,
                               int steps, double x0) {
    /* Track equilibrium from a_start to a_end using delay convention.
     * The jump occurs when the current equilibrium sheet vanishes
     * (i.e., we cross the bifurcation set from inside the cusp to outside).
     *
     * Strategy: follow the nearest equilibrium at each step.
     * When the nearest equilibrium disappears (distance > threshold),
     * we have crossed the fold. */
    double x_cur = x0;
    double a_jump = a_end; /* default: no jump found */

    for (int i = 0; i <= steps; i++) {
        double a = a_start + (a_end - a_start) * i / steps;
        CuspEquilibrium* eq = cusp_equilibrium(a, b_fixed);
        if (!eq || eq->n_roots == 0) {
            a_jump = a;
            cusp_equilibrium_free(eq);
            break;
        }

        /* Find nearest root to x_cur */
        int best = 0;
        double min_dist = fabs(eq->x[0] - x_cur);
        for (int j = 1; j < eq->n_roots; j++) {
            double dist = fabs(eq->x[j] - x_cur);
            if (dist < min_dist) { min_dist = dist; best = j; }
        }

        /* If nearest is much farther than step size, a jump occurred */
        double expected_step = fabs(a_end - a_start) / steps;
        if (min_dist > 10.0 * expected_step + 0.1) {
            a_jump = a;
            cusp_equilibrium_free(eq);
            break;
        }

        x_cur = eq->x[best];
        cusp_equilibrium_free(eq);
    }
    return a_jump;
}

double zeeman_find_maxwell_jump(double a_start, double a_end, double b_fixed,
                                 int steps) {
    /* Find where the global minimum switches between two sheets.
     * For the cusp with b = 0 and sweeping a:
     * Two minima at x = +/- sqrt(-a) for a < 0.
     * They have equal depth (both are global minima).
     *
     * For b != 0: the two minima have different depths.
     * The Maxwell set is b = 0 for a < 0.
     * The Maxwell jump occurs when crossing b = 0. */

    /* For a fixed b, the Maxwell jump in a occurs when the
     * two stable equilibria have equal potential.
     * For the cusp V = x^4/4 + a*x^2/2 + b*x:
     * The minima satisfy x^3 + a*x + b = 0.
     *
     * Equal depth condition: V(x1) = V(x2)
     * Simplifies to: x1^2 + x1*x2 + x2^2 + a = 0 (since x1+x2+x3=0) */

    for (int i = 0; i <= steps; i++) {
        double a = a_start + (a_end - a_start) * i / steps;
        CuspEquilibrium* eq = cusp_equilibrium(a, b_fixed);

        if (eq && eq->n_stable >= 2) {
            /* Find two stable equilibria */
            double xs[2];
            int ns = 0;
            for (int j = 0; j < eq->n_roots && ns < 2; j++) {
                if (eq->stability[j] == STABILITY_STABLE) {
                    xs[ns++] = eq->x[j];
                }
            }
            if (ns >= 2) {
                double V0 = cusp_potential(xs[0], a, b_fixed);
                double V1 = cusp_potential(xs[1], a, b_fixed);
                if ((V0 - V1) * (i > 0 ? 1.0 : -1.0) < 0.0) {
                    /* Sign changed: Maxwell jump occurred */
                    cusp_equilibrium_free(eq);
                    return a;
                }
            }
        }
        cusp_equilibrium_free(eq);
    }
    return a_end; /* no jump found */
}

bool zeeman_detect_bifurcation_crossing(const double* a, const double* b,
                                         int n, int* crossing_idx) {
    for (int i = 1; i < n; i++) {
        double d1 = cusp_discriminant(a[i-1], b[i-1]);
        double d2 = cusp_discriminant(a[i], b[i]);
        if (d1 * d2 < 0.0) {
            if (crossing_idx) *crossing_idx = i;
            return true;
        }
        /* Also check if we hit exactly zero */
        if (fabs(d2) < ZM_EPSILON && fabs(d1) > ZM_EPSILON) {
            if (crossing_idx) *crossing_idx = i;
            return true;
        }
    }
    return false;
}

void zeeman_equilibrium_manifold_section(double a_min, double a_max,
                                          double b_fixed, int n,
                                          double* x_stable_upper,
                                          double* x_stable_lower,
                                          double* x_unstable,
                                          int* n_su, int* n_sl, int* n_unst) {
    *n_su = 0; *n_sl = 0; *n_unst = 0;

    for (int i = 0; i < n; i++) {
        double a = a_min + (a_max - a_min) * i / (n - 1);
        CuspEquilibrium* eq = cusp_equilibrium(a, b_fixed);
        if (!eq) continue;

        for (int j = 0; j < eq->n_roots; j++) {
            if (eq->stability[j] == STABILITY_STABLE) {
                if (eq->x[j] > 0.0 && x_stable_upper) {
                    x_stable_upper[(*n_su)++] = eq->x[j];
                } else if (x_stable_lower) {
                    x_stable_lower[(*n_sl)++] = eq->x[j];
                }
            } else if (eq->stability[j] == STABILITY_UNSTABLE && x_unstable) {
                x_unstable[(*n_unst)++] = eq->x[j];
            }
        }
        cusp_equilibrium_free(eq);
    }
}

/* ==============================================================
 * Divergence Test (Zeeman, 1976)
 *
 * One of the five qualities of catastrophe: small changes in
 * initial conditions lead to large differences in final states.
 *
 * We test divergence by sampling initial conditions along a line
 * in state space and measuring the spread of final states after
 * evolving through the cusp dynamics.
 * ============================================================== */

DivergenceTest* zeeman_divergence_test(double a_param, double b_start,
                                        double b_end, int n_trajectories) {
    DivergenceTest* dt = calloc(1, sizeof(DivergenceTest));
    if (!dt) return NULL;
    dt->n_trajectories = n_trajectories;
    dt->x_final = malloc((size_t)n_trajectories * sizeof(double));
    if (!dt->x_final) { free(dt); return NULL; }

    for (int i = 0; i < n_trajectories; i++) {
        double b = b_start + (b_end - b_start) * i / (n_trajectories - 1);
        double x = 0.0;
        /* Evolve: dx/dt = -(x^3 + a*x + b) */
        double dt_sim = 0.01;
        int steps = 500;
        for (int s = 0; s < steps; s++) {
            double dx = -(x * x * x + a_param * x + b);
            x += dx * dt_sim;
            if (x > 5.0) x = 5.0;
            if (x < -5.0) x = -5.0;
        }
        dt->x_final[i] = x;
    }

    /* Compute divergence measure: variance of final states */
    double mean = 0.0;
    for (int i = 0; i < n_trajectories; i++) mean += dt->x_final[i];
    mean /= n_trajectories;

    double variance = 0.0;
    for (int i = 0; i < n_trajectories; i++) {
        double diff = dt->x_final[i] - mean;
        variance += diff * diff;
    }
    variance /= n_trajectories;
    dt->divergence_measure = sqrt(variance);

    /* Sensitivity: max separation / parameter spread */
    double x_min = dt->x_final[0], x_max = dt->x_final[0];
    for (int i = 1; i < n_trajectories; i++) {
        if (dt->x_final[i] < x_min) x_min = dt->x_final[i];
        if (dt->x_final[i] > x_max) x_max = dt->x_final[i];
    }
    double param_spread = fabs(b_end - b_start);
    dt->sensitivity = (param_spread > ZM_EPSILON)
                    ? (x_max - x_min) / param_spread : 0.0;

    return dt;
}

void zeeman_divergence_test_free(DivergenceTest* dt) {
    if (!dt) return;
    free(dt->x_final);
    free(dt);
}

void zeeman_divergence_test_print(const DivergenceTest* dt) {
    printf("DivergenceTest: n=%d  measure=%.4f  sensitivity=%.4f\n",
           dt->n_trajectories, dt->divergence_measure, dt->sensitivity);
    printf("  Final states range: [%.4f, %.4f]\n",
           dt->x_final[0], dt->x_final[dt->n_trajectories-1]);
}

/* ==============================================================
 * Inaccessibility & Classification
 * ============================================================== */

bool zeeman_is_inaccessible(double a, double b, double x) {
    /* A state x is inaccessible if it lies on the unstable equilibrium
     * sheet (3x^2 + a < 0) of the cusp. */
    CuspEquilibrium* eq = cusp_equilibrium(a, b);
    if (!eq) return false;

    bool result = false;
    for (int i = 0; i < eq->n_roots; i++) {
        if (fabs(eq->x[i] - x) < 0.01
            && eq->stability[i] == STABILITY_UNSTABLE) {
            result = true;
            break;
        }
    }
    cusp_equilibrium_free(eq);
    return result;
}

void zeeman_map_to_cusp_normal(double c3, double c2, double c1, double c0,
                                double* a_out, double* b_out, double* shift) {
    /* Map a general cubic potential V(x) = c3*x^4/4 + c2*x^3/3 + c1*x^2/2 + c0*x
     * to cusp normal form V(y) = y^4/4 + a*y^2/2 + b*y
     * via translation and scaling.
     *
     * Step 1: Eliminate cubic term via x = y - c2/(3*c3)
     * Step 2: Scale y -> y/cbrt(c3) to get coefficient 1/4 on y^4
     *
     * After translation x = y + d:
     *   c3*(y+d)^4/4 = c3*(y^4 + 4dy^3 + ...)/4
     *
     * Need d = -c2/(3*c3) to eliminate cubic term. */
    if (fabs(c3) < ZM_EPSILON) {
        /* Degenerate: not a quartic, can't map to cusp */
        *a_out = 0.0;
        *b_out = 0.0;
        *shift = 0.0;
        return;
    }

    double d = -c2 / (3.0 * c3);
    *shift = d;

    /* Compute a, b in normal form:
     * V(x) expanded around d:
     * V(d+y) = c3*y^4/4 + (c3*3d^2 + c2*d + c1)*y^2/2 + (higher)*y
     *
     * After scaling y -> y / cbrt(c3):
     * a = (3*c3*d^2 + 2*c2*d + c1) / c3^(2/3)
     * b = (c3*d^3 + c2*d^2 + c1*d + c0) / c3^(1/3) */

    double s = pow(c3, 1.0/3.0);
    double s2 = s * s;

    /* Coefficient of y^2 term after translation */
    double coeff_y2 = 3.0 * c3 * d * d + 2.0 * c2 * d + c1;
    *a_out = coeff_y2 / s2;

    /* Coefficient of y term after translation */
    double coeff_y = c3 * d * d * d + c2 * d * d + c1 * d + c0;
    *b_out = coeff_y / s;
}

/* ==============================================================
 * Potential Surface
 * ============================================================== */

double* zeeman_potential_surface(double a_min, double a_max, int na,
                                  double b_min, double b_max, int nb,
                                  double x) {
    double* surface = malloc((size_t)(na * nb) * sizeof(double));
    if (!surface) return NULL;

    double da = (na > 1) ? (a_max - a_min) / (na - 1) : 0.0;
    double db = (nb > 1) ? (b_max - b_min) / (nb - 1) : 0.0;

    for (int ia = 0; ia < na; ia++) {
        double a = a_min + ia * da;
        for (int ib = 0; ib < nb; ib++) {
            double b = b_min + ib * db;
            surface[ia * nb + ib] = cusp_potential(x, a, b);
        }
    }
    return surface;
}

double zeeman_global_minimum(double a, double b) {
    CuspEquilibrium* eq = cusp_equilibrium(a, b);
    if (!eq || eq->n_roots == 0) return 0.0;

    double x_best = eq->x[0];
    double V_best = cusp_potential(x_best, a, b);
    for (int i = 1; i < eq->n_roots; i++) {
        double Vi = cusp_potential(eq->x[i], a, b);
        if (Vi < V_best) { V_best = Vi; x_best = eq->x[i]; }
    }
    cusp_equilibrium_free(eq);
    return x_best;
}

double zeeman_barrier_height(double a, double b) {
    /* Energy barrier between two stable equilibria */
    CuspEquilibrium* eq = cusp_equilibrium(a, b);
    if (!eq || eq->n_stable < 2) {
        cusp_equilibrium_free(eq);
        return 0.0;
    }

    /* Find the two stable equilibria and the unstable one between them */
    double xs[2], V_s[2];
    int ns = 0;
    double xu = 0.0, V_u = 0.0;
    bool has_unstable = false;

    for (int i = 0; i < eq->n_roots; i++) {
        if (eq->stability[i] == STABILITY_STABLE && ns < 2) {
            xs[ns] = eq->x[i];
            V_s[ns] = cusp_potential(xs[ns], a, b);
            ns++;
        } else if (eq->stability[i] == STABILITY_UNSTABLE) {
            xu = eq->x[i];
            V_u = cusp_potential(xu, a, b);
            has_unstable = true;
        }
    }
    cusp_equilibrium_free(eq);

    if (ns < 2 || !has_unstable) return 0.0;

    /* Barrier = V(unstable) - max(V(stable1), V(stable2)) */
    double V_ref = (V_s[0] > V_s[1]) ? V_s[0] : V_s[1];
    return V_u - V_ref;
}

/* ==============================================================
 * Slow-Fast Cusp Dynamics
 *
 * epsilon * dx/dt = -(x^3 + a(x,y)*x + b(x,y))   (fast)
 * dy/dt = g(x, y)                                   (slow)
 *
 * where a(x,y) = a0 + alpha*y
 *       b(x,y) = b0
 *       g(x,y) = beta*(x - gamma*y)
 *
 * This produces relaxation oscillations when the slow variable
 * drives the fast variable back and forth across the cusp folds.
 * ============================================================== */

SlowFastCusp* slowfast_cusp_create(double epsilon, double a0, double b0,
                                    double alpha, double beta, double gamma) {
    SlowFastCusp* sf = calloc(1, sizeof(SlowFastCusp));
    if (!sf) return NULL;
    sf->x = 0.5;
    sf->y = 0.0;
    sf->epsilon = epsilon;
    sf->a0 = a0;
    sf->b0 = b0;
    sf->alpha = alpha;
    sf->beta = beta;
    sf->gamma = gamma;
    sf->t = 0.0;
    return sf;
}

void slowfast_cusp_free(SlowFastCusp* sf) {
    free(sf);
}

void slowfast_cusp_rhs(const SlowFastCusp* sf, double* dx, double* dy) {
    double a = sf->a0 + sf->alpha * sf->y;
    double b = sf->b0;

    /* Fast dynamics */
    *dx = -(sf->x * sf->x * sf->x + a * sf->x + b) / sf->epsilon;

    /* Slow dynamics: y relaxes toward x/gamma */
    *dy = sf->beta * (sf->x - sf->gamma * sf->y);
}

void slowfast_cusp_step(SlowFastCusp* sf, double dt) {
    double dx, dy;
    slowfast_cusp_rhs(sf, &dx, &dy);
    sf->x += dx * dt;
    sf->y += dy * dt;
    sf->t += dt;

    if (sf->x > 3.0)  sf->x = 3.0;
    if (sf->x < -3.0) sf->x = -3.0;
    if (sf->y > 2.0)  sf->y = 2.0;
    if (sf->y < -2.0) sf->y = -2.0;
}

void slowfast_cusp_simulate(SlowFastCusp* sf, double duration, double dt) {
    int steps = (int)(duration / dt);
    for (int i = 0; i < steps; i++) {
        slowfast_cusp_step(sf, dt);
    }
}

bool slowfast_cusp_is_relaxation(const SlowFastCusp* sf) {
    /* Relaxation oscillations occur when epsilon << 1 and
     * the slow variable drives the system across cusp folds */
    return sf->epsilon < 0.2
        && sf->a0 < 0.0  /* bistable when y=0 */
        && fabs(sf->alpha) > 0.5;  /* significant coupling */
}

double slowfast_cusp_period(SlowFastCusp* sf, double duration, double dt) {
    int steps = (int)(duration / dt);
    int crossings = 0;
    double last_cross_time = 0.0;
    bool was_negative = (sf->x < 0.0);

    SlowFastCusp* sf2 = slowfast_cusp_create(sf->epsilon, sf->a0, sf->b0,
                                              sf->alpha, sf->beta, sf->gamma);
    if (!sf2) return 0.0;
    sf2->x = sf->x; sf2->y = sf->y; sf2->t = sf->t;

    for (int i = 0; i < steps; i++) {
        slowfast_cusp_step(sf2, dt);
        bool is_negative = (sf2->x < 0.0);
        if (is_negative != was_negative) {
            double interval = sf2->t - last_cross_time;
            if (crossings > 0 && interval > 0.01) {
                /* crossing detected */;
            }
            last_cross_time = sf2->t;
            crossings++;
        }
        was_negative = is_negative;
    }
    double period = (crossings >= 2) ? (sf2->t / (crossings/2.0)) : 0.0;
    slowfast_cusp_free(sf2);
    return period;
}

void slowfast_cusp_print(const SlowFastCusp* sf) {
    printf("SlowFastCusp: x=%.4f y=%.4f eps=%.4f a0=%.2f b0=%.2f alpha=%.2f t=%.2f\n",
           sf->x, sf->y, sf->epsilon, sf->a0, sf->b0, sf->alpha, sf->t);
    double a = sf->a0 + sf->alpha * sf->y;
    printf("  Effective a=%.4f  Relaxation: %s  Bistable: %s\n",
           a,
           slowfast_cusp_is_relaxation(sf) ? "YES" : "NO",
           cusp_is_bistable(a, sf->b0) ? "YES" : "NO");
}

/* ==============================================================
 * Zeeman's Five Qualities of Catastrophe
 *
 * 1. Bimodality: system has two distinct stable states
 * 2. Sudden jumps: discontinuous change in state
 * 3. Hysteresis: path-dependent behavior
 * 4. Divergence: small cause, large effect
 * 5. Inaccessibility: intermediate states are unstable
 * ============================================================== */

bool zeeman_quality_bimodality(double a, double b) {
    CuspEquilibrium* eq = cusp_equilibrium(a, b);
    bool result = (eq && eq->n_stable >= 2);
    cusp_equilibrium_free(eq);
    return result;
}

bool zeeman_quality_sudden_jump(const double* a, const double* x, int n,
                                 double* jump_a, double* jump_magnitude) {
    if (n < 2) return false;
    double max_dx = 0.0;
    int jump_idx = 0;
    for (int i = 1; i < n; i++) {
        double da = fabs(a[i] - a[i-1]);
        double dx = fabs(x[i] - x[i-1]);
        /* Normalize by parameter step */
        double normalized_dx = (da > ZM_EPSILON) ? dx / da : dx;
        if (normalized_dx > max_dx) {
            max_dx = normalized_dx;
            jump_idx = i;
        }
    }
    if (max_dx > 5.0) {  /* significant normalized jump */
        if (jump_a) *jump_a = a[jump_idx];
        if (jump_magnitude) *jump_magnitude = fabs(x[jump_idx] - x[jump_idx-1]);
        return true;
    }
    return false;
}

double zeeman_quality_hysteresis(double a_min, double a_max, double b, int steps) {
    HysteresisLoop* loop = zeeman_hysteresis_loop(a_min, a_max, b, steps,
                                                    CONVENTION_DELAY);
    if (!loop) return 0.0;
    double area = zeeman_hysteresis_area(loop);
    zeeman_hysteresis_free(loop);
    return area;
}

double zeeman_quality_divergence(double a, double b_start, double b_end, int n) {
    DivergenceTest* dt = zeeman_divergence_test(a, b_start, b_end, n);
    if (!dt) return 0.0;
    double div = dt->divergence_measure;
    zeeman_divergence_test_free(dt);
    return div;
}

bool zeeman_quality_inaccessibility(double a, double b) {
    /* Check if there exists an unstable equilibrium */
    CuspEquilibrium* eq = cusp_equilibrium(a, b);
    if (!eq) return false;
    for (int i = 0; i < eq->n_roots; i++) {
        if (eq->stability[i] == STABILITY_UNSTABLE) {
            cusp_equilibrium_free(eq);
            return true;
        }
    }
    cusp_equilibrium_free(eq);
    return false;
}

void zeeman_five_qualities_report(double a, double b) {
    printf("========================================\n");
    printf(" Zeeman's Five Qualities of Catastrophe\n");
    printf(" Parameters: a=%.4f, b=%.4f\n", a, b);
    printf("========================================\n");

    /* 1. Bimodality */
    bool bimodal = zeeman_quality_bimodality(a, b);
    printf(" 1. Bimodality:        %s", bimodal ? "YES" : "NO");
    if (bimodal) {
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        int count = 0;
        for (int i = 0; i < eq->n_roots; i++) {
            if (eq->stability[i] == STABILITY_STABLE && count < 2) {
                printf(" (stable: x=%.3f)", eq->x[i]);
                count++;
            }
        }
        cusp_equilibrium_free(eq);
    }
    printf("\n");

    /* 2. Sudden jumps */
    int n_test = 100;
    double* a_arr = malloc((size_t)n_test * sizeof(double));
    double* x_arr = malloc((size_t)n_test * sizeof(double));
    for (int i = 0; i < n_test; i++) {
        a_arr[i] = a + 2.0 * ((double)i / (n_test-1) - 0.5);
        CuspEquilibrium* eq = cusp_equilibrium(a_arr[i], b);
        x_arr[i] = (eq && eq->n_roots > 0) ? eq->x[0] : 0.0;
        cusp_equilibrium_free(eq);
    }
    double jump_a, jump_mag;
    bool has_jump = zeeman_quality_sudden_jump(a_arr, x_arr, n_test, &jump_a, &jump_mag);
    printf(" 2. Sudden jumps:      %s", has_jump ? "YES" : "NO");
    if (has_jump) printf(" (at a=%.3f, magnitude=%.3f)", jump_a, jump_mag);
    printf("\n");
    free(a_arr); free(x_arr);

    /* 3. Hysteresis */
    double hyst = zeeman_quality_hysteresis(a - 1.0, a + 1.0, b, 50);
    printf(" 3. Hysteresis:        area=%.4f %s\n",
           hyst, hyst > 0.01 ? "(significant)" : "(negligible)");

    /* 4. Divergence */
    double div = zeeman_quality_divergence(a, b - 0.5, b + 0.5, 20);
    printf(" 4. Divergence:        measure=%.4f %s\n",
           div, div > 0.1 ? "(sensitive)" : "(insensitive)");

    /* 5. Inaccessibility */
    bool inacc = zeeman_quality_inaccessibility(a, b);
    printf(" 5. Inaccessibility:   %s\n", inacc ? "YES (unstable branch exists)" : "NO");

    printf("========================================\n");
}