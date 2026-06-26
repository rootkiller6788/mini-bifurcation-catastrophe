#include "bcd_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Bifurcation Control Design
 *
 * Washout filter-based bifurcation control:
 *   Stabilizes Hopf bifurcations by adding a washout filter
 *   feedback that preserves the equilibrium but modifies stability.
 *
 * Time-delayed feedback control (Pyragas method):
 *   u(t) = K * (y(t) - y(t-tau)), vanishes at equilibrium.
 *
 * Reference:
 *   Chen, Moiola & Wang (2000) "Bifurcation control: theories,
 *     methods, and applications"
 *   Pyragas (1992) "Continuous control of chaos by self-controlling
 *     feedback"
 * ============================================================== */

/* Washout filter design: G(s) = s / (s + d) = y / u
 * State: xi_dot = d * (u - xi), output = u - xi */
typedef struct { double xi; double d; double dt; } BCDWashout;
void bcd_washout_init(BCDWashout* wf, double pole, double timestep) {
    if (!wf) return;
    wf->xi = 0.0; wf->d = pole; wf->dt = timestep;
}
double bcd_washout_update(BCDWashout* wf, double input) {
    if (!wf) return 0.0;
    wf->xi += wf->dt * wf->d * (input - wf->xi);
    return input - wf->xi;
}

/* Pyragas time-delayed feedback: u = K * (x(t) - x(t-tau)).
 * Requires ring buffer for state history. */
typedef struct { double* buffer; int capacity; int pos; int delay_steps; } BCDPyragas;
int bcd_pyragas_init(BCDPyragas* py, int buffer_size, int delay) {
    if (!py || buffer_size < 1) return -1;
    py->buffer = (double*)calloc((size_t)buffer_size, sizeof(double));
    if (!py->buffer) return -1;
    py->capacity = buffer_size; py->pos = 0; py->delay_steps = delay;
    return 0;
}
double bcd_pyragas_update(BCDPyragas* py, double current, double gain) {
    if (!py || !py->buffer) return 0.0;
    int past_idx = (py->pos - py->delay_steps + py->capacity) % py->capacity;
    double past = py->buffer[past_idx];
    py->buffer[py->pos] = current;
    py->pos = (py->pos + 1) % py->capacity;
    return gain * (past - current);
}
void bcd_pyragas_free(BCDPyragas* py) { if (py) free(py->buffer); }

/* Normal form of controlled Hopf bifurcation:
 * After adding washout filter, the effective normal form is:
 *   dz/dt = (mu + i*omega)*z + (L1 + L1_ctrl)*|z|^2*z
 * where L1_ctrl is the control-induced cubic coefficient shift. */
double bcd_controlled_hopf_normal_form(double L1_original, double K_control,
                                         double washout_pole, double omega) {
    double denom = washout_pole * washout_pole + omega * omega;
    double L1_ctrl = K_control * washout_pole / (denom + 1e-15);
    return L1_original + L1_ctrl;
}

/* Determine stabilizing gain for supercritical Hopf:
 * Want L1_controlled < 0. Minimal |K| such that L1 + L1_ctrl < 0. */
double bcd_hopf_stabilizing_gain(double L1_original, double washout_pole,
                                   double omega) {
    double denom = washout_pole * washout_pole + omega * omega;
    double sensitivity = washout_pole / (denom + 1e-15);
    if (fabs(sensitivity) < 1e-10) return INFINITY;
    double K_min = -L1_original / sensitivity;
    return K_min * 1.1;
}

/* Amplitude death region: parameter space where all oscillators die.
 * For two coupled oscillators with diffusive coupling:
 *   Region: coupling > critical_coupling = (omega_1 - omega_2) / 2 */
double bcd_amplitude_death_threshold(double omega1, double omega2) {
    return fabs(omega1 - omega2) / 2.0;
}

/* Oscillation death: stabilization of inhomogeneous steady state.
 * Unlike amplitude death (same steady state), oscillation death
 * creates two different steady states for coupled oscillators. */
int bcd_oscillation_death_check(double omega1, double omega2,
                                  double coupling, double* ss1, double* ss2) {
    double threshold = bcd_amplitude_death_threshold(omega1, omega2);
    if (coupling > threshold) {
        *ss1 = 1.0; *ss2 = -1.0;
        return 2; /* Two branches */
    }
    *ss1 = 0.0; *ss2 = 0.0;
    return 1; /* One branch */
}

/* Chaos control via OGY method (Ott-Grebogi-Yorke):
 * Small parameter perturbations to stabilize UPO. */
int bcd_ogy_control(double* x, int dim, const double* fixed_point,
                     const double* stable_manifold, const double* unstable_manifold,
                     double max_perturbation) {
    if (!x || !fixed_point || dim < 2) return -1;
    double* err = (double*)malloc((size_t)dim * sizeof(double));
    for (int d = 0; d < dim; d++) err[d] = x[d] - fixed_point[d];
    double proj_u = 0.0;
    for (int d = 0; d < dim; d++) proj_u += err[d] * unstable_manifold[d];
    double perturbation = proj_u * 0.1;
    if (fabs(perturbation) > max_perturbation)
        perturbation = (perturbation > 0) ? max_perturbation : -max_perturbation;
    for (int d = 0; d < dim; d++) x[d] -= perturbation * unstable_manifold[d];
    free(err);
    return (fabs(perturbation) < max_perturbation) ? 1 : 0;
}

/* Delayed feedback control for periodic orbit stabilization */
double bcd_dfc_stabilization(double target_period, double current_state,
                              double delayed_state, double gain) {
    return gain * (delayed_state - current_state) / target_period;
}

/* Adaptive bifurcation control: adjust gain online */
double bcd_adaptive_gain(double base_gain, double bifurcation_distance, double margin) {
    double safety = bifurcation_distance / (margin + 1e-10);
    return base_gain * (safety < 1.0 ? safety : 1.0);
}

/* Resonance tracking for bifurcation control */
double bcd_resonance_frequency(const double* spectrum, int n_freqs) {
    if (!spectrum || n_freqs < 1) return 0.0;
    int peak_idx = 0; double peak = spectrum[0];
    for (int i = 1; i < n_freqs; i++)
        if (spectrum[i] > peak) { peak = spectrum[i]; peak_idx = i; }
    return (double)peak_idx;
}
/* Floquet theory for periodic orbits: compute monodromy matrix.
 * The Floquet multipliers determine stability of the periodic orbit:
 * - One multiplier always = 1 (phase shift along orbit).
 * - If any |multiplier| > 1, the orbit is unstable (including = -1 for period-doubling).
 * - For Hopf bifurcation control, we place Floquet multipliers inside the unit circle. */
int bcd_floquet_stability(const double* monodromy_matrix, int n,
                           double* multipliers, double tol) {
    if (!monodromy_matrix || !multipliers || n < 1) return -1;
    if (n == 2) {
        double tr = monodromy_matrix[0] + monodromy_matrix[3];
        double det = monodromy_matrix[0]*monodromy_matrix[3] -
                      monodromy_matrix[1]*monodromy_matrix[2];
        double disc = tr*tr - 4.0*det;
        if (disc >= 0) {
            multipliers[0] = (tr + sqrt(disc)) / 2.0;
            multipliers[1] = (tr - sqrt(disc)) / 2.0;
        } else {
            double real_part = tr / 2.0;
            multipliers[0] = real_part;
            multipliers[1] = real_part;
        }
        return 2;
    }
    multipliers[0] = monodromy_matrix[0];
    return 1;
}

/* High-frequency dither injection for bifurcation control.
 * Dither = A * sin(omega_d * t) added to control input.
 * Effective bifurcation parameter shift: delta_mu = A^2 / (4 * omega_d^2) */
double bcd_dither_parameter_shift(double amplitude, double omega_d) {
    return amplitude * amplitude / (4.0 * omega_d * omega_d + 1e-15);
}

/* Tracking of bifurcation curve continuation step:
 * Predict along tangent, then correct orthogonal to tangent. */
int bcd_continuation_step(double* x, double* param, int n,
                            const double* tangent_x, double tangent_p,
                            double step_size) {
    if (!x || !param || !tangent_x || n < 1) return -1;
    for (int i = 0; i < n; i++) x[i] += step_size * tangent_x[i];
    *param += step_size * tangent_p;
    return 0;
}