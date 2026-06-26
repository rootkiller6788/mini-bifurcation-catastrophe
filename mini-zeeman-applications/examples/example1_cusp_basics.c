#include "zeeman_core.h"
#include "zeeman_dynamics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=======================================\n");
    printf("  Zeeman Cusp Catastrophe — Demo 1\n");
    printf("  Core Geometry & Hysteresis\n");
    printf("=======================================\n\n");

    /* 1. Visualize the cusp potential */
    printf("--- 1. Cusp Potential V(x) = x^4/4 + a*x^2/2 + b*x ---\n");
    double b_vals[] = {-0.5, 0.0, 0.5};
    double a_test = -1.0;

    for (int ib = 0; ib < 3; ib++) {
        double b = b_vals[ib];
        printf("\na=%.1f, b=%.1f:\n", a_test, b);
        CuspEquilibrium* eq = cusp_equilibrium(a_test, b);
        cusp_equilibrium_print(eq);
        for (int i = 0; i < eq->n_roots; i++) {
            printf("  V(x=%.3f) = %.6f\n",
                   eq->x[i], cusp_potential(eq->x[i], a_test, b));
        }
        printf("  Bistable: %s\n", cusp_is_bistable(a_test, b) ? "YES" : "NO");
        cusp_equilibrium_free(eq);
    }

    /* 2. Hysteresis loop: sweep a from -2 to 0 at b=0 */
    printf("\n--- 2. Hysteresis Loop (b=0, a: -2 -> 0 -> -2) ---\n");
    HysteresisLoop* loop = zeeman_hysteresis_loop(-2.0, 0.0, 0.0, 40,
                                                   CONVENTION_DELAY);
    zeeman_hysteresis_print(loop);
    printf("\nForward sweep (a=-2 to 0):\n");
    for (int i = 0; i <= 40; i += 8) {
        printf("  a=%.2f -> x=%.4f\n", loop->a_fwd[i], loop->x_fwd[i]);
    }
    printf("Backward sweep (a=0 to -2):\n");
    for (int i = 0; i <= 40; i += 8) {
        printf("  a=%.2f -> x=%.4f\n", loop->a_bwd[i], loop->x_bwd[i]);
    }
    zeeman_hysteresis_free(loop);

    /* 3. Bifurcation diagram */
    printf("\n--- 3. Bifurcation Diagram ---\n");
    BifurcationDiagram* bd = zeeman_bifurcation_diagram(-2.0, 1.0, -0.5, 0.5, 30, 30);
    zeeman_bifurcation_diagram_print(bd);
    printf("Sample of equilibrium manifold:\n");
    for (int i = 0; i < 20 && i < bd->n_points; i += 5) {
        printf("  (a=%.3f, b=%.3f) x=%.4f V=%.4f stab=%d\n",
               bd->points[i].a, bd->points[i].b,
               bd->points[i].x, bd->points[i].V,
               (int)bd->points[i].stability);
    }
    zeeman_bifurcation_diagram_free(bd);

    /* 4. Zeeman machine simulation */
    printf("\n--- 4. Zeeman Catastrophe Machine ---\n");
    ZeemanMachine* zm = zeeman_machine_create_default();
    printf("Initial state:\n");
    zeeman_machine_print(zm);

    printf("\nSweeping control peg C from (0,0) to (-2,2):\n");
    ZeemanPath* path = zeeman_path_create(50, CONVENTION_DELAY);
    zeeman_machine_sweep_control(zm, 0.0, 0.0, -2.0, 2.0, 40, path);
    zeeman_path_print(path);
    int j_idx;
    if (zeeman_path_detect_jump(path, &j_idx)) {
        printf("  JUMP detected at index %d: theta %.3f -> %.3f\n",
               j_idx, path->x[j_idx-1], path->x[j_idx]);
    }
    zeeman_path_free(path);
    zeeman_machine_free(zm);

    /* 5. Five qualities */
    printf("\n--- 5. Zeeman Five Qualities ---\n");
    zeeman_five_qualities_report(-1.0, 0.0);

    printf("\n=======================================\n");
    printf("  Demo complete.\n");
    printf("=======================================\n");
    return 0;
}