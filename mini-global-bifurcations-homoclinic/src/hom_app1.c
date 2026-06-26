/* hom_app1.c -- L7: Ship capsize homoclinic tangles (Thompson 1989).
 * The escape equation: xddot + beta*xdot + x - x^2 = F*sin(omega*t)
 * models a ship rolling in beam seas. The unperturbed system has a
 * homoclinic orbit connecting the saddle at the potential barrier.
 * Melnikov analysis predicts capsize threshold vs wave amplitude. */
#include "../include/homoclinic_melnikov.h"
#include "../include/homoclinic_core.h"
#include <math.h>
#include <stdio.h>

static void ship_unperturbed(HOMState x, HOMState* dx, void* p) {
    (void)p;
    dx->x[0] = x.x[1];  /* angle_dot */
    dx->x[1] = -x.x[0] + x.x[0] * x.x[0];  /* xddot = -x + x^2 */
}
static double wave_forcing(double t, void* p) {
    double* pars = (double*)p;
    return pars[0] * sin(pars[1] * t);  /* F*sin(omega*t) */
}

int hom_ship_capsize(double F_wave, double omega, double damping,
                     double* capsize_threshold, double* safe_heading) {
    double pars[] = {F_wave, omega};
    HOMMelnikov* mel = hom_melnikov_create(ship_unperturbed, NULL, wave_forcing, pars);

    /* Generate unperturbed homoclinic orbit (analytic approx) */
    int n_orbit = 300;
    HOMState* orbit = (HOMState*)calloc((size_t)n_orbit, sizeof(HOMState));
    double t_start = -6.0, t_end = 6.0, dt = (t_end - t_start) / (double)(n_orbit - 1);
    for (int i = 0; i < n_orbit; i++) {
        double t = t_start + i * dt;
        double sech2 = 1.0 / (cosh(t) * cosh(t));
        orbit[i].x[0] = 1.5 * sech2;
        orbit[i].x[1] = -3.0 * sech2 * tanh(t);
        orbit[i].dim = 2;
    }

    double M0 = hom_melnikov_integral(mel, orbit, n_orbit, 0.0, dt);
    double thresh = hom_melnikov_bifurcation_threshold(mel);
    *capsize_threshold = thresh;

    printf("[hom_app1] Ship capsize: F=%.3f omega=%.2f damping=%.3f\n", F_wave, omega, damping);
    printf("  Melnikov M(t0=0)=%.6f thr=%.6f\n", M0, thresh);
    printf("  Capsize risk: %s\n", F_wave > thresh ? "HIGH (homoclinic tangles)" : "LOW");

    bool chaos = hom_melnikov_chaos_prediction(mel, damping);
    printf("  Chaos prediction: %s\n", chaos ? "CHAOTIC rolling" : "Regular motion");

    *safe_heading = (thresh - F_wave) / F_wave * 100.0;

    hom_melnikov_free(mel);
    free(orbit);
    return 0;
}

int main(void) {
    double thr, safe;
    hom_ship_capsize(0.8, 1.5, 0.1, &thr, &safe);
    printf("=== hom_app1: Ship Capsize Homoclinic ===\n");
    printf("threshold=%.4f safe_margin=%.1f%%\n", thr, safe);
    return 0;
}
