#ifndef ZEEMAN_CORE_H
#define ZEEMAN_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* Cusp catastrophe geometry, Zeeman machine, core types */

#define ZM_PI            3.14159265358979323846
#define ZM_E             2.71828182845904523536
#define ZM_EPSILON       1e-12
#define ZM_GOLDEN_RATIO  1.6180339887498948482
#define ZM_SQRT2         1.4142135623730950488
#define ZM_SQRT3         1.7320508075688772935
#define ZM_DEG2RAD       (ZM_PI / 180.0)
#define ZM_RAD2DEG       (180.0 / ZM_PI)
#define ZM_INV_PI        0.31830988618379067154
#define ZM_TWO_PI        6.28318530717958647692
#define ZM_HALF_PI       1.57079632679489661923

#define ZM_DEFAULT_RADIUS     1.0
#define ZM_DEFAULT_K1         1.0
#define ZM_DEFAULT_K2         1.0
#define ZM_DEFAULT_FIXED_X    2.0
#define ZM_DEFAULT_FIXED_Y    0.0
#define ZM_DEFAULT_DAMPING    0.5

typedef enum {
    CATASTROPHE_NONE        = 0,
    CATASTROPHE_FOLD        = 1,
    CATASTROPHE_CUSP        = 2,
    CATASTROPHE_SWALLOWTAIL = 3,
    CATASTROPHE_BUTTERFLY   = 4,
    CATASTROPHE_HYPERBOLIC  = 5,
    CATASTROPHE_ELLIPTIC    = 6,
    CATASTROPHE_PARABOLIC   = 7
} CatastropheType;

typedef enum {
    STABILITY_UNKNOWN    = 0,
    STABILITY_STABLE     = 1,
    STABILITY_UNSTABLE   = 2,
    STABILITY_DEGENERATE = 3,
    STABILITY_METASTABLE = 4
} StabilityType;

typedef enum {
    CONVENTION_DELAY      = 0,
    CONVENTION_MAXWELL    = 1,
    CONVENTION_STOCHASTIC = 2
} JumpConvention;

typedef struct { double a; double b; } CuspParams;

typedef struct {
    double x[3];
    StabilityType stability[3];
    int    n_roots;
    int    n_stable;
    double discriminant;
} CuspEquilibrium;

typedef struct {
    double R, k1, k2;
    double Fx, Fy, Cx, Cy;
    double theta, theta_dot, damping, inertia;
} ZeemanMachine;

typedef struct { double x, a, b, V; StabilityType stability; } CuspSurfacePoint;

typedef struct {
    double *a, *b, *x;
    int n, cap;
    JumpConvention convention;
} ZeemanPath;

typedef struct {
    double *a_fwd, *x_fwd, *a_bwd, *x_bwd;
    int n;
    double jump_fwd_a, jump_bwd_a, hysteresis_area;
} HysteresisLoop;

typedef struct {
    CuspSurfacePoint *points;
    int n_points, cap;
    double a_min, a_max, b_min, b_max;
} BifurcationDiagram;

/* Cusp catastrophe core math */
double  cusp_potential(double x, double a, double b);
double  cusp_force(double x, double a, double b);
double  cusp_curvature(double x, double a);
CuspEquilibrium* cusp_equilibrium(double a, double b);
void             cusp_equilibrium_free(CuspEquilibrium* eq);
double           cusp_discriminant(double a, double b);
bool             cusp_is_bistable(double a, double b);
bool             cusp_is_on_bifurcation_set(double a, double b);
void             cusp_bifurcation_set_point(double t, double* a, double* b);
double           cusp_maxwell_b(double a);
int              cusp_num_equilibria(double a, double b);
void             cusp_equilibrium_print(const CuspEquilibrium* eq);

/* Higher catastrophes */
double cusp_potential_fold(double x, double a);
double cusp_potential_swallowtail(double x, double a, double b, double c);
double cusp_potential_butterfly(double x, double a, double b, double c, double d);

/* Zeeman machine */
ZeemanMachine* zeeman_machine_create(double R, double k1, double k2, double Fx, double Fy);
ZeemanMachine* zeeman_machine_create_default(void);
void           zeeman_machine_free(ZeemanMachine* m);
void           zeeman_machine_set_control(ZeemanMachine* m, double Cx, double Cy);
void           zeeman_machine_get_control(const ZeemanMachine* m, double* Cx, double* Cy);
double         zeeman_machine_get_theta(const ZeemanMachine* m);
void           zeeman_machine_set_theta(ZeemanMachine* m, double theta);
void           zeeman_machine_set_damping(ZeemanMachine* m, double damping);
double         zeeman_machine_potential(const ZeemanMachine* m, double theta);
double         zeeman_machine_torque(const ZeemanMachine* m, double theta);
CuspEquilibrium* zeeman_machine_equilibrium(const ZeemanMachine* m);
void           zeeman_machine_to_cusp(const ZeemanMachine* m, double* a, double* b);
void           zeeman_machine_from_cusp(ZeemanMachine* m, double a, double b, double scale);
void           zeeman_machine_step(ZeemanMachine* m, double dt);
void           zeeman_machine_sweep_control(ZeemanMachine* m, double Cx_start, double Cy_start, double Cx_end, double Cy_end, int steps, ZeemanPath* path_out);
void           zeeman_machine_print(const ZeemanMachine* m);
double         zeeman_machine_band_length(const ZeemanMachine* m, double peg_x, double peg_y, double theta);

/* Path tracking */
ZeemanPath* zeeman_path_create(int capacity, JumpConvention convention);
void        zeeman_path_free(ZeemanPath* path);
void        zeeman_path_add(ZeemanPath* path, double a, double b, double x);
void        zeeman_path_track_equilibrium(ZeemanPath* path, double a_start, double b_start, double a_end, double b_end, int steps);
void        zeeman_path_print(const ZeemanPath* path);
bool        zeeman_path_detect_jump(const ZeemanPath* path, int* jump_idx);
double      zeeman_path_jump_magnitude(const ZeemanPath* path);
ZeemanPath* zeeman_plan_safe_path(double a_start, double b_start, double a_end, double b_end, int steps);

/* Hysteresis */
HysteresisLoop* zeeman_hysteresis_loop(double a_start, double a_end, double b_fixed, int steps, JumpConvention convention);
void            zeeman_hysteresis_free(HysteresisLoop* loop);
double          zeeman_hysteresis_area(const HysteresisLoop* loop);
void            zeeman_hysteresis_print(const HysteresisLoop* loop);

/* Bifurcation diagram */
BifurcationDiagram* zeeman_bifurcation_diagram(double a_min, double a_max, double b_min, double b_max, int na, int nb);
void                zeeman_bifurcation_diagram_free(BifurcationDiagram* bd);
int                 zeeman_bifurcation_diagram_count_catastrophes(const BifurcationDiagram* bd);
void                zeeman_bifurcation_diagram_print(const BifurcationDiagram* bd);

/* Numerical utilities */
int  zeeman_solve_cubic_cardano(double a, double b, double* roots);
double zeeman_cubic_real_root(double a, double b);
CatastropheType zeeman_classify_catastrophe(double a, double b, double c, double d);

#endif /* ZEEMAN_CORE_H */