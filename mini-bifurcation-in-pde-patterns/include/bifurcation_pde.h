#ifndef BIFURCATION_PDE_H
#define BIFURCATION_PDE_H
#include "pattern_formation.h"
typedef struct{double param;double*eigenvalues;int n_ev;}BifurcationPoint;
typedef struct{double* params;double**branches;int n_params;int n_states;}BifurcationDiagram;
BifurcationPoint* bp_create(int n_ev);
void bp_free(BifurcationPoint* bp);
BifurcationDiagram* bd_create(int n_params,int n_states);
void bd_free(BifurcationDiagram* bd);
void bd_add_point(BifurcationDiagram* bd,int idx,double param,double* state);
void bp_linear_stability_1d(const RDField* f,TuringParams tp,int* n_unstable);
void bp_continuation_step(RDField* f,TuringParams* tp,double d_param);
int bp_detect_bifurcation(const double* eigenvalues,int n,double tol);
void bp_compute_nullspace(const double* A,int N,double* null_vec);
#define BP_MAX_EV 64
#define BP_TOL 1e-6
#endif
