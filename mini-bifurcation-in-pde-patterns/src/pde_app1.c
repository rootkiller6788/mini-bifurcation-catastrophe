/* pde_app1.c -- L7: Turing pattern in reaction-diffusion (Turing 1952).
 * The Schnakenberg model: du/dt = Du*d2u/dx2 + gamma*(a - u + u^2*v)
 *                        dv/dt = Dv*d2v/dx2 + gamma*(b - u^2*v)
 * Spontaneous pattern formation when Turing instability condition
 * is met. Compute dispersion relation and dominant wavelength. */
#include "../include/reaction_diffusion.h"
#include "../include/numerical_pde.h"
#include <math.h>
#include <stdio.h>

int pde_zebra_pattern(double Du, double Dv, double a, double b,
                      double gamma, double* wavelength, double* growth_rate) {
    /* Steady state: u0 = a + b, v0 = b/(a+b)^2 */
    double u0 = a + b;
    double v0 = b / (u0 * u0);
    double fu = gamma * (-1.0 + 2.0 * u0 * v0);
    double fv = gamma * u0 * u0;
    double gu = gamma * (-2.0 * u0 * v0);
    double gv = gamma * (-u0 * u0);

    printf("[pde_app1] Turing pattern (Schnakenberg): Du=%.2f Dv=%.2f a=%.2f b=%.2f\n",
           Du, Dv, a, b);
    printf("  Steady state: u0=%.4f v0=%.4f\n", u0, v0);
    printf("  Jacobian: fu=%.4f fv=%.4f gu=%.4f gv=%.4f\n", fu, fv, gu, gv);

    /* Turing condition: fu + gv < 0, fu*gv - fv*gu > 0, Dv*fu + Du*gv > 2*sqrt(Du*Dv*(fu*gv-fv*gu)) */
    double det = fu * gv - fv * gu;
    double trace = fu + gv;
    double turing_cond = Dv * fu + Du * gv;
    double turing_th = 2.0 * sqrt(Du * Dv * det);

    printf("  Turing conditions:\n");
    printf("    fu+gv=%.4f < 0: %s\n", trace, trace < 0 ? "YES" : "NO");
    printf("    det=%.4f > 0: %s\n", det, det > 0 ? "YES" : "NO");
    printf("    Dv*fu+Du*gv=%.4f > 2*sqrt(...)=%.4f: %s\n",
           turing_cond, turing_th, turing_cond > turing_th ? "YES (Turing!)" : "NO");

    /* Dominant wavenumber: k_max^2 = sqrt(det/(Du*Dv)) */
    double k2_max = sqrt(det / (Du * Dv));
    *wavelength = 2.0 * M_PI / sqrt(k2_max);
    *growth_rate = det - (k2_max * k2_max * Du * Dv);

    printf("  Dominant wavelength: %.4f, growth rate: %.6f\n", *wavelength, *growth_rate);

    /* 1D simulation on grid */
    int N = 128;
    Grid1D* g = (Grid1D*)calloc(1, sizeof(Grid1D));
    g->N = N; g->dx = *wavelength * 4.0 / (double)N;

    /* Dispersion relation */
    printf("  Dispersion (k, sigma):\n");
    for (int i = 0; i <= 10; i++) {
        double k = 2.0 * M_PI * i / (*wavelength * 4.0);
        double sigma = 0.5 * (fu + gv - k*k*(Du + Dv) +
                        sqrt((fu + gv - k*k*(Du+Dv))*(fu+gv - k*k*(Du+Dv)) -
                             4.0*(det - k*k*(Dv*fu + Du*gv) + k*k*k*k*Du*Dv)));
        if (i % 2 == 0) printf("    k=%.4f sigma=%.6f\n", k, sigma);
    }

    free(g);
    return 0;
}

int main(void) {
    double wl, gr;
    pde_zebra_pattern(0.05, 1.0, 0.1, 0.9, 100.0, &wl, &gr);
    printf("=== pde_app1: Turing Pattern ===\n");
    printf("wavelength=%.4f growth=%.6f\n", wl, gr);
    return 0;
}
