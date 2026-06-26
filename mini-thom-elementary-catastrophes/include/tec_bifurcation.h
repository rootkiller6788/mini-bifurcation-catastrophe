#ifndef TEC_BIFURCATION_H
#define TEC_BIFURCATION_H
#include "tec_core.h"

/* ==============================================================
 * tec_bifurcation.h - Bifurcation Set & Catastrophe Manifold
 *
 * The bifurcation set is the projection of the singularity set
 * onto the control space. It divides control space into regions
 * with different numbers of equilibria.
 *
 * For fold: single curve a = 0
 * For cusp: semi-cubical parabola 4a^3 + 27b^2 = 0
 * For swallowtail: surface defined by quartic discriminant
 * For butterfly: hypersurface defined by quintic discriminant
 *
 * Maxwell set: locus where two or more minima have equal depth
 * (relevant for systems with hysteresis).
 * ============================================================== */

/* --- Bifurcation Set --- */
typedef struct {
    TEC_BifurcationPoint* points;
    int n_points;
    int capacity;
    TEC_CatastropheType type;
} TEC_BifurcationSet;

/* --- Bifurcation Diagram Slice --- */
typedef struct {
    double* param_range;     /* Varying parameter values */
    double* equilibria;      /* Equilibrium positions (may have multiple per param) */
    int* n_eq_per_param;
    int n_params;
    int max_eq;
} TEC_BifurcationDiagram;

/* --- API: Bifurcation Set Computation --- */
TEC_BifurcationSet* tec_bif_set_fold(int n_pts);
TEC_BifurcationSet* tec_bif_set_cusp(int n_pts, double a_range);
TEC_BifurcationSet* tec_bif_set_swallowtail(int n_pts);
void tec_bif_set_free(TEC_BifurcationSet* bs);
void tec_bif_set_print(const TEC_BifurcationSet* bs);
bool tec_is_on_bifurcation_set(TEC_CatastropheType type, const double* params);

/* --- API: Bifurcation Diagram --- */
TEC_BifurcationDiagram* tec_bif_diagram_fold(int n_pts, double a_range);
TEC_BifurcationDiagram* tec_bif_diagram_cusp(int n_pts, double b_range, double a_fixed);
void tec_bif_diagram_free(TEC_BifurcationDiagram* bd);
void tec_bif_diagram_print(const TEC_BifurcationDiagram* bd);

/* --- API: Maxwell Set --- */
TEC_MaxwellPoint* tec_maxwell_set_cusp(int n_pts, double a_range);
TEC_MaxwellPoint* tec_maxwell_set_butterfly(double a, double b, double c);
void tec_maxwell_free(TEC_MaxwellPoint* mp, int n);
void tec_maxwell_print(const TEC_MaxwellPoint* mp);

/* --- API: Hysteresis Analysis --- */
bool tec_has_hysteresis(TEC_CatastropheType type, const double* params);
double tec_hysteresis_loop_width(TEC_CatastropheType type, const double* params);
int tec_count_equilibria_region(TEC_CatastropheType type, const double* params);

/* --- API: Delay Convention --- */
double tec_delay_convention_jump(TEC_CatastropheType type, const double* params, double x_start, double param_step);
double tec_bif_distance_to_set(TEC_CatastropheType t, const double* params);
int tec_bif_crossing_detect(TEC_CatastropheType t, const double* p1, const double* p2);

#endif
