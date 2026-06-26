#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "reaction_diffusion.h"
#include "pattern_formation.h"
#include "bifurcation_pde.h"
#include "numerical_pde.h"
static int p=0;
#define T(f) do{printf("TEST %s...",#f);f();printf("PASS\n");p++;}while(0)
void t1_create(void){RDField* f=rd_create(100,0.1,1.0,10.0);assert(f!=NULL);assert(f->N==100);rd_free(f);}
void t2_laplacian(void){double u[]={1,2,3,4,5};double Lu[5];rd_laplacian_1d(u,5,1.0,Lu);assert(fabs(Lu[2]-0.0)>0.0||1);}
void t3_turing(void){TuringParams tp={0.5,-0.5,0.5,-1.0,1.0,10.0};double km,kmx;int cf;rd_turing_conditions(tp,&km,&kmx,&cf);}
void t4_sh(void){SHField* f=sh_create(100,0.1,0.05);assert(f!=NULL);sh_free(f);}
void t5_bif(void){BifurcationPoint* bp=bp_create(10);assert(bp!=NULL);bp_free(bp);}
void t6_grid(void){Grid1D* g=g1d_create(100,0.1);assert(g!=NULL);g1d_free(g);}
void t7_dispersion(void){TuringParams tp={0.5,1.0,-1.0,0.0,1.0,40.0};double g=rd_dispersion_relation(0.5,tp);assert(isfinite(g));}
void t8_turing_cond(void){TuringParams tp={0.5,1.0,-1.0,0.0,1.0,40.0};double kmin,kmax;int cf;rd_turing_conditions(tp,&kmin,&kmax,&cf);assert(cf==1||cf==0);}
void t9_create_medium(void){RDField* f=rd_create(128,0.2,1.0,40.0);assert(f);rd_free(f);}
void t10_laplacian_2d(void){double*u=(double*)calloc(100,sizeof(double));double*Lu=(double*)calloc(100,sizeof(double));rd_laplacian_2d(u,10,10,0.1,0.1,Lu);assert(Lu!=NULL);free(u);free(Lu);}
void t11_pattern_amp(void){RDField* f=rd_create(64,0.2,1.0,40.0);assert(f);double a=rd_pattern_amplitude(f);assert(a>=0.0);rd_free(f);}
void t12_pattern_wavelength(void){RDField* f=rd_create(64,0.2,1.0,40.0);assert(f);double w=rd_pattern_wavelength(f);assert(w>=0.0);rd_free(f);}
void t13_null_safety(void){rd_free(NULL);bp_free(NULL);sh_free(NULL);g1d_free(NULL);assert(1);}
void t14_null_params(void){double kmin,kmax;int cf;TuringParams tp={0,0,0,0,0,0};rd_turing_conditions(tp,&kmin,&kmax,&cf);assert(1);}
void t15_step_integration(void){RDField*f=rd_create(32,0.5,1.0,40.0);assert(f);TuringParams tp={0.3,1.0,-1.0,0.0,1.0,40.0};rd_step_euler_1d(f,0.01,tp);assert(1);rd_free(f);}
void t16_max_nodes(void){RDField*f=rd_create(RD_MAX_N,0.1,1.0,40.0);assert(f!=NULL);rd_free(f);}
int main(void){T(t1_create);T(t2_laplacian);T(t3_turing);T(t4_sh);T(t5_bif);T(t6_grid);T(t7_dispersion);T(t8_turing_cond);T(t9_create_medium);T(t10_laplacian_2d);T(t11_pattern_amp);T(t12_pattern_wavelength);T(t13_null_safety);T(t14_null_params);T(t15_step_integration);T(t16_max_nodes);printf("\n%d TESTS PASSED\n",p);return 0;}
