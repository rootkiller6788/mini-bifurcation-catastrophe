#ifndef PATTERN_FORMATION_H
#define PATTERN_FORMATION_H
#include "reaction_diffusion.h"
typedef struct{double*field;int N;double dx,epsilon;}SHField;
typedef struct{double A_re,A_im;double omega;}AmplitudeMode;
typedef struct{AmplitudeMode*modes;int n_modes;double kc;}AmplitudeEq;
SHField* sh_create(int N,double dx,double epsilon);
void sh_free(SHField* f);
void sh_step_euler(SHField* f,double dt);
void sh_step_etdrk4(SHField* f,double dt);
void sh_linear_stability(double k,double epsilon,double*sigma);
double sh_most_unstable_wavenumber(double epsilon);
void sh_compute_pattern(SHField* f,double*amplitude,double*wavelength);
AmplitudeEq* ae_create(int n_modes,double kc);
void ae_free(AmplitudeEq* ae);
void ae_step(AmplitudeEq* ae,double dt);
double ae_landau_coefficient(double kc,double epsilon);
#define SH_DEFAULT_EPS 0.1
#define SH_MAX_N 1024
#endif
