#include "zeeman_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Cusp Catastrophe - Core Mathematics
 *
 * V(x; a, b) = x^4/4 + a*x^2/2 + b*x
 *
 * This is the universal unfolding of the singularity x^4.
 * Every 4-parameter family of potentials with a degenerate
 * critical point can be reduced to this form (Thom, 1975).
 * ============================================================== */

double cusp_potential(double x, double a, double b) {
    double x2 = x * x;
    return 0.25 * x2 * x2 + 0.5 * a * x2 + b * x;
}

double cusp_force(double x, double a, double b) {
    /* -dV/dx = -(x^3 + a*x + b) */
    return -(x * x * x + a * x + b);
}

double cusp_curvature(double x, double a) {
    /* d^2V/dx^2 = 3*x^2 + a */
    return 3.0 * x * x + a;
}

/* ==============================================================
 * Discriminant and Bifurcation Set
 *
 * For the depressed cubic x^3 + a*x + b = 0:
 *   discriminant delta = b^2/4 + a^3/27
 *
 *   delta > 0: 1 real root
 *   delta < 0: 3 distinct real roots
 *   delta = 0: multiple root (bifurcation set)
 *
 * Bifurcation set equation: 4*a^3 + 27*b^2 = 0
 * Equivalent to delta = 0 (multiplying by 108).
 * ============================================================== */

double cusp_discriminant(double a, double b) {
    return b * b / 4.0 + a * a * a / 27.0;
}

bool cusp_is_bistable(double a, double b) {
    /* Bistable when inside the cusp: a < 0 and delta < 0 */
    if (a >= 0.0) return false;
    return cusp_discriminant(a, b) < 0.0;
}

bool cusp_is_on_bifurcation_set(double a, double b) {
    double d = 4.0 * a * a * a + 27.0 * b * b;
    return fabs(d) < ZM_EPSILON;
}

void cusp_bifurcation_set_point(double t, double* a, double* b) {
    /* Parametric form: a = -3*t^2, b = 2*t^3 for t in R */
    *a = -3.0 * t * t;
    *b = 2.0 * t * t * t;
}

double cusp_maxwell_b(double a) {
    /* Maxwell set: b = 0 for a < 0 (two minima equal depth) */
    (void)a;
    return 0.0;
}

int cusp_num_equilibria(double a, double b) {
    double d = cusp_discriminant(a, b);
    if (d > ZM_EPSILON) return 1;
    if (d < -ZM_EPSILON) return 3;
    /* On bifurcation set: check if triple root */
    if (fabs(a) < ZM_EPSILON && fabs(b) < ZM_EPSILON) return 1; /* triple at origin */
    return 2; /* one simple, one double */
}

/* ==============================================================
 * Cardano's Method for the Depressed Cubic x^3 + a*x + b = 0
 *
 * Step 1: Compute discriminant delta = b^2/4 + a^3/27
 *
 * Step 2: If delta > 0 (one real root):
 *   u = cbrt(-b/2 + sqrt(delta))
 *   v = cbrt(-b/2 - sqrt(delta))
 *   x1 = u + v (real)
 *
 * Step 3: If delta < 0 (three real roots, casus irreducibilis):
 *   r = sqrt(-a^3/27) = sqrt(-a/3)^3
 *   phi = acos(-b/(2*r))
 *   x1 = 2 * cbrt(r) * cos(phi/3)
 *   x2 = 2 * cbrt(r) * cos((phi + 2*pi)/3)
 *   x3 = 2 * cbrt(r) * cos((phi + 4*pi)/3)
 *
 * Step 4: If delta = 0 (multiple root):
 *   x1 = 2 * cbrt(-b/2)
 *   x2 = x3 = -cbrt(-b/2)
 * ============================================================== */

int zeeman_solve_cubic_cardano(double a, double b, double* roots) {
    double delta = cusp_discriminant(a, b);
    int n;

    if (delta > ZM_EPSILON) {
        /* One real root */
        double sqrt_delta = sqrt(delta);
        double u = -b / 2.0 + sqrt_delta;
        double v = -b / 2.0 - sqrt_delta;
        double cu = (u >= 0.0) ? pow(u, 1.0/3.0) : -pow(-u, 1.0/3.0);
        double cv = (v >= 0.0) ? pow(v, 1.0/3.0) : -pow(-v, 1.0/3.0);
        roots[0] = cu + cv;
        n = 1;
    } else if (delta < -ZM_EPSILON) {
        /* Three real roots (casus irreducibilis) */
        double r = sqrt(-a * a * a / 27.0);
        double phi = acos(-b / (2.0 * r));
        double cr = 2.0 * pow(r, 1.0/3.0);
        roots[0] = cr * cos(phi / 3.0);
        roots[1] = cr * cos((phi + 2.0 * ZM_PI) / 3.0);
        roots[2] = cr * cos((phi + 4.0 * ZM_PI) / 3.0);
        n = 3;
    } else {
        /* Multiple root on bifurcation set */
        if (fabs(b) < ZM_EPSILON && fabs(a) < ZM_EPSILON) {
            roots[0] = 0.0;
            n = 1; /* triple root */
        } else {
            double cb = (b >= 0.0) ? pow(b/2.0, 1.0/3.0) : -pow(-b/2.0, 1.0/3.0);
            roots[0] = 2.0 * cb;
            roots[1] = -cb;
            n = 2;
        }
    }
    return n;
}

double zeeman_cubic_real_root(double a, double b) {
    double roots[3];
    zeeman_solve_cubic_cardano(a, b, roots);
    return roots[0];
}

/* ==============================================================
 * CuspEquilibrium - Full equilibrium analysis
 * ============================================================== */

CuspEquilibrium* cusp_equilibrium(double a, double b) {
    CuspEquilibrium* eq = calloc(1, sizeof(CuspEquilibrium));
    if (!eq) return NULL;

    eq->discriminant = cusp_discriminant(a, b);
    double roots[3];
    eq->n_roots = zeeman_solve_cubic_cardano(a, b, roots);
    eq->n_stable = 0;

    for (int i = 0; i < eq->n_roots; i++) {
        eq->x[i] = roots[i];
        double curv = cusp_curvature(roots[i], a);
        if (curv > ZM_EPSILON) {
            eq->stability[i] = STABILITY_STABLE;
            eq->n_stable++;
        } else if (curv < -ZM_EPSILON) {
            eq->stability[i] = STABILITY_UNSTABLE;
        } else {
            eq->stability[i] = STABILITY_DEGENERATE;
        }
    }

    /* Sort: stable first, then unstable, then degenerate */
    for (int i = 0; i < eq->n_roots - 1; i++) {
        for (int j = i + 1; j < eq->n_roots; j++) {
            if ((int)eq->stability[i] > (int)eq->stability[j]) {
                double tx = eq->x[i];
                eq->x[i] = eq->x[j];
                eq->x[j] = tx;
                StabilityType ts = eq->stability[i];
                eq->stability[i] = eq->stability[j];
                eq->stability[j] = ts;
            }
        }
    }

    return eq;
}

void cusp_equilibrium_free(CuspEquilibrium* eq) {
    free(eq);
}

void cusp_equilibrium_print(const CuspEquilibrium* eq) {
    printf("CuspEquilibrium: a=... n_roots=%d n_stable=%d delta=%.6e\n",
           eq->n_roots, eq->n_stable, eq->discriminant);
    for (int i = 0; i < eq->n_roots; i++) {
        const char* stab_str;
        switch (eq->stability[i]) {
            case STABILITY_STABLE:     stab_str = "STABLE"; break;
            case STABILITY_UNSTABLE:   stab_str = "UNSTABLE"; break;
            case STABILITY_DEGENERATE: stab_str = "DEGENERATE"; break;
            case STABILITY_METASTABLE: stab_str = "METASTABLE"; break;
            default:                   stab_str = "UNKNOWN"; break;
        }
        printf("  root[%d]: x=% .6f  V=%.6f  %s\n",
               i, eq->x[i], cusp_potential(eq->x[i], 0.0, 0.0), stab_str);
    }
}

/* ==============================================================
 * Higher-Order Catastrophes
 * ============================================================== */

double cusp_potential_fold(double x, double a) {
    /* V(x; a) = x^3/3 + a*x */
    return x * x * x / 3.0 + a * x;
}

double cusp_potential_swallowtail(double x, double a, double b, double c) {
    /* V(x; a, b, c) = x^5/5 + a*x^3/3 + b*x^2/2 + c*x */
    double x2 = x * x;
    return x2 * x2 * x / 5.0 + a * x2 * x / 3.0 + b * x2 / 2.0 + c * x;
}

double cusp_potential_butterfly(double x, double a, double b, double c, double d) {
    /* V(x; a,b,c,d) = x^6/6 + a*x^4/4 + b*x^3/3 + c*x^2/2 + d*x */
    double x2 = x * x;
    return x2 * x2 * x2 / 6.0 + a * x2 * x2 / 4.0
         + b * x2 * x / 3.0 + c * x2 / 2.0 + d * x;
}

CatastropheType zeeman_classify_catastrophe(double a, double b, double c, double d) {
    /* Classify by codimension and control parameter count */
    (void)a; (void)b; (void)c; (void)d;
    /* Simplest classification: check which parameters are effectively non-zero */
    double thresh = ZM_EPSILON;
    if (fabs(c) < thresh && fabs(d) < thresh) {
        if (fabs(b) < thresh) {
            return (fabs(a) < thresh) ? CATASTROPHE_NONE : CATASTROPHE_FOLD;
        }
        return CATASTROPHE_CUSP;
    }
    if (fabs(d) < thresh) return CATASTROPHE_SWALLOWTAIL;
    return CATASTROPHE_BUTTERFLY;
}

/* ==============================================================
 * Zeeman Catastrophe Machine
 *
 * Physical model: disk of radius R at origin, free to rotate.
 * Pin P(theta) = R*(cos theta, sin theta) on disk perimeter.
 * Elastic band 1: fixed peg F to P.  Energy = 0.5*k1*|P-F|^2
 * Elastic band 2: control peg C to P.  Energy = 0.5*k2*|P-C|^2
 *
 * Total: V(theta) = 0.5*k1*((R*cos-Fx)^2+(R*sin-Fy)^2)
 *                 + 0.5*k2*((R*cos-Cx)^2+(R*sin-Cy)^2)
 *
 * Expanded: V(theta) = const_term + alpha*cos(theta) + beta*sin(theta)
 *
 * where:
 *   alpha = -R*(k1*Fx + k2*Cx)
 *   beta  = -R*(k1*Fy + k2*Cy)
 *
 * Minima occur where tan(theta) = beta/alpha.
 *
 * As C moves, the equilibrium theta changes continuously
 * until a fold catastrophe occurs (two equilibria merge).
 * ============================================================== */

ZeemanMachine* zeeman_machine_create(double R, double k1, double k2,
                                      double Fx, double Fy) {
    ZeemanMachine* m = calloc(1, sizeof(ZeemanMachine));
    if (!m) return NULL;
    m->R  = R;
    m->k1 = k1;
    m->k2 = k2;
    m->Fx = Fx;
    m->Fy = Fy;
    m->Cx = 0.0;
    m->Cy = 0.0;
    m->theta = 0.0;
    m->theta_dot = 0.0;
    m->damping = ZM_DEFAULT_DAMPING;
    m->inertia = 1.0;
    return m;
}

ZeemanMachine* zeeman_machine_create_default(void) {
    return zeeman_machine_create(ZM_DEFAULT_RADIUS, ZM_DEFAULT_K1, ZM_DEFAULT_K2,
                                  ZM_DEFAULT_FIXED_X, ZM_DEFAULT_FIXED_Y);
}

void zeeman_machine_free(ZeemanMachine* m) {
    free(m);
}

void zeeman_machine_set_control(ZeemanMachine* m, double Cx, double Cy) {
    m->Cx = Cx;
    m->Cy = Cy;
}

void zeeman_machine_get_control(const ZeemanMachine* m, double* Cx, double* Cy) {
    *Cx = m->Cx;
    *Cy = m->Cy;
}

double zeeman_machine_get_theta(const ZeemanMachine* m) {
    return m->theta;
}

void zeeman_machine_set_theta(ZeemanMachine* m, double theta) {
    m->theta = theta;
}

void zeeman_machine_set_damping(ZeemanMachine* m, double damping) {
    m->damping = damping;
}

/* ==============================================================
 * Machine Potential and Torque
 *
 * V(theta) = 0.5*k1*|P - F|^2 + 0.5*k2*|P - C|^2
 *
 * |P - F|^2 = (R*cos - Fx)^2 + (R*sin - Fy)^2
 *           = R^2 + Fx^2 + Fy^2 - 2*R*(Fx*cos + Fy*sin)
 *
 * So V(theta) = const + alpha*cos(theta) + beta*sin(theta)
 *
 * Torque = -dV/dtheta = alpha*sin(theta) - beta*cos(theta)
 * ============================================================== */

double zeeman_machine_band_length(const ZeemanMachine* m,
                                   double peg_x, double peg_y, double theta) {
    double px = m->R * cos(theta);
    double py = m->R * sin(theta);
    double dx = px - peg_x;
    double dy = py - peg_y;
    return sqrt(dx * dx + dy * dy);
}

double zeeman_machine_potential(const ZeemanMachine* m, double theta) {
    double px = m->R * cos(theta);
    double py = m->R * sin(theta);
    double d1x = px - m->Fx, d1y = py - m->Fy;
    double d2x = px - m->Cx, d2y = py - m->Cy;
    double d1_sq = d1x * d1x + d1y * d1y;
    double d2_sq = d2x * d2x + d2y * d2y;
    return 0.5 * m->k1 * d1_sq + 0.5 * m->k2 * d2_sq;
}

double zeeman_machine_torque(const ZeemanMachine* m, double theta) {
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    double alpha = -m->R * (m->k1 * m->Fx + m->k2 * m->Cx);
    double beta  = -m->R * (m->k1 * m->Fy + m->k2 * m->Cy);
    /* torque = -dV/dtheta = -(alpha*(-sin) + beta*cos) = alpha*sin - beta*cos */
    return alpha * sin_t - beta * cos_t;
}

/* ==============================================================
 * Machine Equilibrium: find theta where torque = 0
 *
 * tan(theta) = beta/alpha when alpha != 0
 *
 * Actually the potential V(theta) = A - B*cos(theta - phi)
 * where B = sqrt(alpha^2 + beta^2), phi = atan2(beta, alpha)
 *
 * So V(theta) = const - B*cos(theta - phi)
 * Min at theta = phi, max at theta = phi + pi.
 *
 * But as C moves, B and phi change, and new equilibria can appear.
 * ============================================================== */

CuspEquilibrium* zeeman_machine_equilibrium(const ZeemanMachine* m) {
    CuspEquilibrium* eq = calloc(1, sizeof(CuspEquilibrium));
    if (!eq) return NULL;

    double alpha = -m->R * (m->k1 * m->Fx + m->k2 * m->Cx);
    double beta  = -m->R * (m->k1 * m->Fy + m->k2 * m->Cy);
    double B = sqrt(alpha * alpha + beta * beta);

    if (B < ZM_EPSILON) {
        /* Degenerate: no preferred direction, all angles equivalent */
        eq->n_roots = 0;
        eq->n_stable = 0;
        eq->discriminant = 0.0;
        return eq;
    }

    /* Two equilibria: minimum and maximum */
    double phi = atan2(beta, alpha);
    eq->n_roots = 2;
    eq->n_stable = 1;
    eq->discriminant = -1.0; /* always 2 equilibria in this model */
    eq->x[0] = phi;                    /* stable minimum */
    eq->stability[0] = STABILITY_STABLE;
    eq->x[1] = phi + ZM_PI;            /* unstable maximum */
    eq->stability[1] = STABILITY_UNSTABLE;

    /* Normalize to [0, 2*pi) */
    for (int i = 0; i < 2; i++) {
        while (eq->x[i] < 0.0) eq->x[i] += ZM_TWO_PI;
        while (eq->x[i] >= ZM_TWO_PI) eq->x[i] -= ZM_TWO_PI;
    }

    return eq;
}

/* ==============================================================
 * Map machine coordinates to cusp normal form (a, b)
 *
 * The mapping from control peg (Cx, Cy) to cusp parameters (a, b)
 * involves a nonlinear coordinate transformation.
 *
 * Near the catastrophe point, the machine potential V(theta)
 * is locally equivalent to the cusp catastrophe.
 *
 * From Zeeman (1972): for the standard machine with F=(2,0),
 * the cusp point in control space is near C=(-1, 0), and
 * the bifurcation set forms a cusp-shaped curve.
 *
 * Approximate mapping (after Taylor expansion):
 *   a ~ -(Cx + 1) + higher order
 *   b ~ Cy
 * ============================================================== */

void zeeman_machine_to_cusp(const ZeemanMachine* m, double* a, double* b) {
    /* Approximate mapping from (Cx, Cy) to cusp (a, b)
     * Based on Taylor expansion of V(theta) around the catastrophe point.
     * For the default machine (F = (2,0), R=1, k1=k2=1),
     * the cusp point is approximately at C = (-1, 0). */
    double cx_shifted = m->Cx + 1.0;  /* shift to cusp point */
    double cy = m->Cy;
    /* Nonlinear transformation to cusp parameters */
    *a = -(cx_shifted - 0.5 * cy * cy);
    *b = cy * (1.0 - 0.3 * cx_shifted);
}

void zeeman_machine_from_cusp(ZeemanMachine* m, double a, double b,
                               double scale) {
    /* Inverse mapping from cusp (a,b) to control (Cx, Cy) */
    double cx_shifted = -a * scale;
    double cy = b * scale;
    m->Cx = cx_shifted - 1.0;  /* unshift */
    m->Cy = cy;
}

/* ==============================================================
 * Machine Dynamics
 *
 * Equation: I * d^2(theta)/dt^2 = torque - damping * d(theta)/dt
 * In state-space form:
 *   d(theta)/dt = theta_dot
 *   d(theta_dot)/dt = (torque - damping*theta_dot) / inertia
 *
 * Euler integration for simplicity.
 * ============================================================== */

void zeeman_machine_step(ZeemanMachine* m, double dt) {
    double tau = zeeman_machine_torque(m, m->theta);
    double alpha = (tau - m->damping * m->theta_dot) / m->inertia;
    m->theta_dot += alpha * dt;
    m->theta += m->theta_dot * dt;
    /* Keep theta in reasonable range */
    while (m->theta < -ZM_PI) m->theta += ZM_TWO_PI;
    while (m->theta > ZM_PI)  m->theta -= ZM_TWO_PI;
}

void zeeman_machine_sweep_control(ZeemanMachine* m,
                                   double Cx_start, double Cy_start,
                                   double Cx_end,   double Cy_end,
                                   int steps, ZeemanPath* path_out) {
    if (steps < 1) return;
    for (int i = 0; i <= steps; i++) {
        double t = (double)i / steps;
        double Cx = Cx_start + t * (Cx_end - Cx_start);
        double Cy = Cy_start + t * (Cy_end - Cy_start);
        zeeman_machine_set_control(m, Cx, Cy);
        /* Let machine settle to equilibrium */
        for (int j = 0; j < 100; j++) {
            zeeman_machine_step(m, 0.01);
        }
        if (path_out) {
            double a, b;
            zeeman_machine_to_cusp(m, &a, &b);
            zeeman_path_add(path_out, a, b, m->theta);
        }
    }
}

void zeeman_machine_print(const ZeemanMachine* m) {
    printf("ZeemanMachine: R=%.2f k1=%.2f k2=%.2f\n", m->R, m->k1, m->k2);
    printf("  Fixed peg:  F=(%.2f, %.2f)\n", m->Fx, m->Fy);
    printf("  Control peg: C=(%.2f, %.2f)\n", m->Cx, m->Cy);
    printf("  State: theta=%.4f rad (%.1f deg)  theta_dot=%.4f\n",
           m->theta, m->theta * ZM_RAD2DEG, m->theta_dot);
    printf("  Potential at current theta: %.6f\n",
           zeeman_machine_potential(m, m->theta));
    double a, b;
    zeeman_machine_to_cusp(m, &a, &b);
    printf("  Cusp params: a=%.4f b=%.4f  bistable=%s\n",
           a, b, cusp_is_bistable(a, b) ? "YES" : "NO");
}

/* ==============================================================
 * Path Tracking
 * ============================================================== */

ZeemanPath* zeeman_path_create(int capacity, JumpConvention convention) {
    ZeemanPath* p = calloc(1, sizeof(ZeemanPath));
    if (!p) return NULL;
    p->cap = (capacity > 0) ? capacity : 256;
    p->a = malloc((size_t)p->cap * sizeof(double));
    p->b = malloc((size_t)p->cap * sizeof(double));
    p->x = malloc((size_t)p->cap * sizeof(double));
    p->n = 0;
    p->convention = convention;
    if (!p->a || !p->b || !p->x) {
        zeeman_path_free(p);
        return NULL;
    }
    return p;
}

void zeeman_path_free(ZeemanPath* path) {
    if (!path) return;
    free(path->a);
    free(path->b);
    free(path->x);
    free(path);
}

void zeeman_path_add(ZeemanPath* path, double a, double b, double x) {
    if (path->n >= path->cap) {
        path->cap *= 2;
        path->a = realloc(path->a, (size_t)path->cap * sizeof(double));
        path->b = realloc(path->b, (size_t)path->cap * sizeof(double));
        path->x = realloc(path->x, (size_t)path->cap * sizeof(double));
    }
    path->a[path->n] = a;
    path->b[path->n] = b;
    path->x[path->n] = x;
    path->n++;
}

void zeeman_path_track_equilibrium(ZeemanPath* path,
                                    double a_start, double b_start,
                                    double a_end,   double b_end,
                                    int steps) {
    if (steps < 1) return;
    /* Start with equilibrium at initial point */
    CuspEquilibrium* eq_start = cusp_equilibrium(a_start, b_start);
    if (!eq_start || eq_start->n_roots == 0) {
        cusp_equilibrium_free(eq_start);
        return;
    }
    double x_current = eq_start->x[0]; /* take first (stable) root */
    cusp_equilibrium_free(eq_start);

    for (int i = 0; i <= steps; i++) {
        double t = (double)i / steps;
        double a = a_start + t * (a_end - a_start);
        double b = b_start + t * (b_end - b_start);

        /* Find nearest root of current cubic */
        CuspEquilibrium* eq = cusp_equilibrium(a, b);
        if (!eq || eq->n_roots == 0) {
            cusp_equilibrium_free(eq);
            continue;
        }

        /* Select equilibrium based on convention */
        double x_new;
        if (path->convention == CONVENTION_MAXWELL) {
            /* Maxwell: choose global minimum */
            x_new = eq->x[0];
            double V_min = cusp_potential(eq->x[0], a, b);
            for (int j = 1; j < eq->n_roots; j++) {
                double Vj = cusp_potential(eq->x[j], a, b);
                if (Vj < V_min) {
                    V_min = Vj;
                    x_new = eq->x[j];
                }
            }
        } else {
            /* Delay convention: stay on nearest sheet */
            x_new = eq->x[0];
            double min_dist = fabs(eq->x[0] - x_current);
            for (int j = 1; j < eq->n_roots; j++) {
                double dist = fabs(eq->x[j] - x_current);
                if (dist < min_dist) {
                    min_dist = dist;
                    x_new = eq->x[j];
                }
            }
        }

        zeeman_path_add(path, a, b, x_new);
        x_current = x_new;
        cusp_equilibrium_free(eq);
    }
}

void zeeman_path_print(const ZeemanPath* path) {
    printf("ZeemanPath: %d points, convention=%d\n", path->n, (int)path->convention);
    if (path->n == 0) return;
    printf("  Start: (a=%.4f, b=%.4f, x=%.4f)\n", path->a[0], path->b[0], path->x[0]);
    printf("  End:   (a=%.4f, b=%.4f, x=%.4f)\n",
           path->a[path->n-1], path->b[path->n-1], path->x[path->n-1]);
}

bool zeeman_path_detect_jump(const ZeemanPath* path, int* jump_idx) {
    if (path->n < 2) return false;
    double max_jump = 0.0;
    int idx = 0;
    for (int i = 1; i < path->n; i++) {
        double dx = fabs(path->x[i] - path->x[i-1]);
        if (dx > max_jump) {
            max_jump = dx;
            idx = i;
        }
    }
    if (max_jump > 0.5) { /* threshold for "jump" */
        if (jump_idx) *jump_idx = idx;
        return true;
    }
    return false;
}

double zeeman_path_jump_magnitude(const ZeemanPath* path) {
    if (path->n < 2) return 0.0;
    double max_jump = 0.0;
    for (int i = 1; i < path->n; i++) {
        double dx = fabs(path->x[i] - path->x[i-1]);
        if (dx > max_jump) max_jump = dx;
    }
    return max_jump;
}

ZeemanPath* zeeman_plan_safe_path(double a_start, double b_start,
                                   double a_end,   double b_end, int steps) {
    /* Plan a path that avoids the bifurcation set
     * by going around the cusp region. */
    ZeemanPath* path = zeeman_path_create(steps + 1, CONVENTION_DELAY);

    /* Check if direct path crosses cusp interior */
    bool cross_cusp = false;
    for (int i = 0; i <= 10; i++) {
        double t = i / 10.0;
        double a = a_start + t * (a_end - a_start);
        double b = b_start + t * (b_end - b_start);
        if (cusp_is_bistable(a, b)) { cross_cusp = true; break; }
    }

    if (!cross_cusp) {
        /* Direct path is safe */
        zeeman_path_track_equilibrium(path, a_start, b_start, a_end, b_end, steps);
    } else {
        /* Detour: go above the cusp region via positive b */
        double a_mid = (a_start + a_end) / 2.0;
        double b_safe = fmax(fabs(b_start), fabs(b_end)) + 2.0;
        zeeman_path_track_equilibrium(path, a_start, b_start, a_mid, b_safe, steps/3);
        zeeman_path_track_equilibrium(path, a_mid, b_safe, a_mid, b_safe, steps/3);
        zeeman_path_track_equilibrium(path, a_mid, b_safe, a_end, b_end, steps/3);
    }
    return path;
}

/* ==============================================================
 * Hysteresis Analysis
 * ============================================================== */

HysteresisLoop* zeeman_hysteresis_loop(double a_start, double a_end,
                                        double b_fixed, int steps,
                                        JumpConvention convention) {
    HysteresisLoop* loop = calloc(1, sizeof(HysteresisLoop));
    if (!loop) return NULL;
    loop->n = steps + 1;

    loop->a_fwd = malloc((size_t)loop->n * sizeof(double));
    loop->x_fwd = malloc((size_t)loop->n * sizeof(double));
    loop->a_bwd = malloc((size_t)loop->n * sizeof(double));
    loop->x_bwd = malloc((size_t)loop->n * sizeof(double));

    if (!loop->a_fwd || !loop->x_fwd || !loop->a_bwd || !loop->x_bwd) {
        zeeman_hysteresis_free(loop);
        return NULL;
    }

    /* Forward sweep: a_start -> a_end */
    CuspEquilibrium* eq_start = cusp_equilibrium(a_start, b_fixed);
    double x_cur = (eq_start && eq_start->n_roots > 0) ? eq_start->x[0] : 0.0;
    cusp_equilibrium_free(eq_start);

    loop->jump_fwd_a = a_end; /* default: no jump */

    for (int i = 0; i <= steps; i++) {
        double a = a_start + (a_end - a_start) * i / steps;
        loop->a_fwd[i] = a;

        CuspEquilibrium* eq = cusp_equilibrium(a, b_fixed);
        if (eq && eq->n_roots > 0) {
            /* Stay on nearest sheet (delay convention) */
            int best = 0;
            double min_dist = fabs(eq->x[0] - x_cur);
            for (int j = 1; j < eq->n_roots; j++) {
                double dist = fabs(eq->x[j] - x_cur);
                if (dist < min_dist) { min_dist = dist; best = j; }
            }
            /* Check if previous sheet vanished (jump) */
            if (min_dist > 1.0 && eq->n_roots >= 1) {
                best = 0;
                loop->jump_fwd_a = a;
            }
            x_cur = eq->x[best];
            loop->x_fwd[i] = x_cur;
        } else {
            loop->x_fwd[i] = x_cur;
        }
        cusp_equilibrium_free(eq);
    }

    /* Backward sweep: a_end -> a_start */
    eq_start = cusp_equilibrium(a_end, b_fixed);
    x_cur = (eq_start && eq_start->n_roots > 0) ? eq_start->x[0] : 0.0;
    cusp_equilibrium_free(eq_start);

    loop->jump_bwd_a = a_start; /* default: no jump */

    for (int i = 0; i <= steps; i++) {
        double a = a_end + (a_start - a_end) * i / steps;
        loop->a_bwd[i] = a;

        CuspEquilibrium* eq = cusp_equilibrium(a, b_fixed);
        if (eq && eq->n_roots > 0) {
            int best = 0;
            double min_dist = fabs(eq->x[0] - x_cur);
            for (int j = 1; j < eq->n_roots; j++) {
                double dist = fabs(eq->x[j] - x_cur);
                if (dist < min_dist) { min_dist = dist; best = j; }
            }
            if (min_dist > 1.0 && eq->n_roots >= 1) {
                best = (eq->n_roots > 1) ? eq->n_roots - 1 : 0;
                loop->jump_bwd_a = a;
            }
            x_cur = eq->x[best];
            loop->x_bwd[i] = x_cur;
        } else {
            loop->x_bwd[i] = x_cur;
        }
        cusp_equilibrium_free(eq);
    }

    /* Compute hysteresis area using trapezoidal rule */
    loop->hysteresis_area = 0.0;
    for (int i = 0; i < steps; i++) {
        double dx = loop->a_fwd[i+1] - loop->a_fwd[i];
        double y_diff = (loop->x_fwd[i] - loop->x_bwd[i]
                       + loop->x_fwd[i+1] - loop->x_bwd[i+1]) / 2.0;
        loop->hysteresis_area += fabs(dx * y_diff);
    }

    return loop;
}

void zeeman_hysteresis_free(HysteresisLoop* loop) {
    if (!loop) return;
    free(loop->a_fwd);
    free(loop->x_fwd);
    free(loop->a_bwd);
    free(loop->x_bwd);
    free(loop);
}

double zeeman_hysteresis_area(const HysteresisLoop* loop) {
    return loop->hysteresis_area;
}

void zeeman_hysteresis_print(const HysteresisLoop* loop) {
    printf("HysteresisLoop: n=%d\n", loop->n);
    printf("  Forward jump at a=%.4f\n", loop->jump_fwd_a);
    printf("  Backward jump at a=%.4f\n", loop->jump_bwd_a);
    printf("  Hysteresis area: %.6f\n", loop->hysteresis_area);
    printf("  Forward:  a[0]=%.3f x[0]=%.3f -> a[n-1]=%.3f x[n-1]=%.3f\n",
           loop->a_fwd[0], loop->x_fwd[0],
           loop->a_fwd[loop->n-1], loop->x_fwd[loop->n-1]);
    printf("  Backward: a[0]=%.3f x[0]=%.3f -> a[n-1]=%.3f x[n-1]=%.3f\n",
           loop->a_bwd[0], loop->x_bwd[0],
           loop->a_bwd[loop->n-1], loop->x_bwd[loop->n-1]);
}

/* ==============================================================
 * Bifurcation Diagram
 * ============================================================== */

BifurcationDiagram* zeeman_bifurcation_diagram(double a_min, double a_max,
                                                double b_min, double b_max,
                                                int na, int nb) {
    BifurcationDiagram* bd = calloc(1, sizeof(BifurcationDiagram));
    if (!bd) return NULL;
    bd->a_min = a_min; bd->a_max = a_max;
    bd->b_min = b_min; bd->b_max = b_max;
    bd->cap = na * nb * 3; /* max 3 equilibria per point */
    bd->points = malloc((size_t)bd->cap * sizeof(CuspSurfacePoint));
    bd->n_points = 0;
    if (!bd->points) { free(bd); return NULL; }

    double da = (na > 1) ? (a_max - a_min) / (na - 1) : 0.0;
    double db = (nb > 1) ? (b_max - b_min) / (nb - 1) : 0.0;

    for (int ia = 0; ia < na; ia++) {
        double a = a_min + ia * da;
        for (int ib = 0; ib < nb; ib++) {
            double b = b_min + ib * db;
            CuspEquilibrium* eq = cusp_equilibrium(a, b);
            if (!eq) continue;
            for (int k = 0; k < eq->n_roots; k++) {
                if (bd->n_points >= bd->cap) {
                    bd->cap *= 2;
                    bd->points = realloc(bd->points,
                        (size_t)bd->cap * sizeof(CuspSurfacePoint));
                }
                CuspSurfacePoint* sp = &bd->points[bd->n_points++];
                sp->x = eq->x[k];
                sp->a = a;
                sp->b = b;
                sp->V = cusp_potential(eq->x[k], a, b);
                sp->stability = eq->stability[k];
            }
            cusp_equilibrium_free(eq);
        }
    }
    return bd;
}

void zeeman_bifurcation_diagram_free(BifurcationDiagram* bd) {
    if (!bd) return;
    free(bd->points);
    free(bd);
}

int zeeman_bifurcation_diagram_count_catastrophes(const BifurcationDiagram* bd) {
    int count = 0;
    for (int i = 0; i < bd->n_points; i++) {
        if (bd->points[i].stability == STABILITY_DEGENERATE) count++;
    }
    return count;
}

void zeeman_bifurcation_diagram_print(const BifurcationDiagram* bd) {
    printf("BifurcationDiagram: %d points\n", bd->n_points);
    printf("  Domain: a=[%.2f, %.2f] b=[%.2f, %.2f]\n",
           bd->a_min, bd->a_max, bd->b_min, bd->b_max);
    int n_stable = 0, n_unstable = 0, n_deg = 0;
    for (int i = 0; i < bd->n_points; i++) {
        switch (bd->points[i].stability) {
            case STABILITY_STABLE:     n_stable++; break;
            case STABILITY_UNSTABLE:   n_unstable++; break;
            case STABILITY_DEGENERATE: n_deg++; break;
            default: break;
        }
    }
    printf("  Stable=%d Unstable=%d Degenerate=%d\n", n_stable, n_unstable, n_deg);
}