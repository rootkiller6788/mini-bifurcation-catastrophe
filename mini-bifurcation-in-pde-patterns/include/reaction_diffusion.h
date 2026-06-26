#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H
#include <stdbool.h>
#include <stddef.h>
typedef struct{double*u,*v;int N;double dx,Du,Dv;}RDField;
typedef struct{double alpha,beta,gamma,delta,Du,Dv;}TuringParams;
typedef struct{double*real,*imag;int N;}ComplexField;
RDField* rd_create(int N,double dx,double Du,double Dv);
void rd_free(RDField* f);
void rd_set_initial(RDField* f,double(*u0)(double),double(*v0)(double));
void rd_laplacian_1d(const double* u,int N,double dx,double* Lu);
void rd_laplacian_2d(const double* u,int Nx,int Ny,double dx,double dy,double* Lu);
void rd_step_euler_1d(RDField* f,double dt,TuringParams tp);
void rd_step_rk4_1d(RDField* f,double dt,TuringParams tp);
void rd_turing_conditions(TuringParams tp,double*k2_min,double*k2_max,int*can_form);
double rd_dispersion_relation(double k2,TuringParams tp);
void rd_compute_fft_1d(const double* u,int N,ComplexField* cf);
void rd_compute_ifft_1d(ComplexField* cf,double* u);
double rd_pattern_wavelength(const RDField* f);
double rd_pattern_amplitude(const RDField* f);
#define RD_MAX_N 2048
#define RD_DEFAULT_DT 0.001
#endif
