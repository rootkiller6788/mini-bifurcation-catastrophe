#ifndef NUMERICAL_PDE_H
#define NUMERICAL_PDE_H
#include "bifurcation_pde.h"
typedef struct{double*data;int N;double dx;}Grid1D;
typedef struct{double*data;int Nx,Ny;double dx,dy;}Grid2D;
Grid1D* g1d_create(int N,double dx);void g1d_free(Grid1D* g);
Grid2D* g2d_create(int Nx,int Ny,double dx,double dy);void g2d_free(Grid2D* g);
void g1d_sin_transform(const double* u,int N,double* u_hat);
void g1d_inv_sin_transform(const double* u_hat,int N,double* u);
void g2d_laplacian_5pt(const double* u,int Nx,int Ny,double dx,double dy,double* Lu);
void g2d_laplacian_9pt(const double* u,int Nx,int Ny,double dx,double dy,double* Lu);
void npde_crank_nicolson_1d(double* u,int N,double dt,double D,double dx);
void npde_implicit_euler_1d(double* u,int N,double dt,double D,double dx);
void npde_etd1_1d(double* u,int N,double dt,double D,double dx);
void npde_etdrk4_1d(double* u,int N,double dt,double D,double dx);
#define NPDE_SOLVER_TOL 1e-8
#endif
