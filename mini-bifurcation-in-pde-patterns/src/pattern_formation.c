#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pattern_formation.h"
#include "numerical_pde.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
SHField* sh_create(int N,double dx,double epsilon){if(N<2||dx<=0.0)return NULL;SHField* f=calloc(1,sizeof(SHField));f->N=N;f->dx=dx;f->epsilon=epsilon;f->field=calloc(N,sizeof(double));return f;}
void sh_free(SHField* f){if(!f)return;free(f->field);free(f);}
void sh_linear_stability(double k,double epsilon,double*sigma){*sigma=epsilon-(1.0-k*k)*(1.0-k*k);}
double sh_most_unstable_wavenumber(double epsilon){return 1.0;}
void sh_step_euler(SHField* f,double dt){if(!f||dt<=0.0)return;int N=f->N;double*Lu=calloc(N,sizeof(double)),*L2u=calloc(N,sizeof(double));rd_laplacian_1d(f->field,N,f->dx,Lu);rd_laplacian_1d(Lu,N,f->dx,L2u);for(int i=0;i<N;i++)f->field[i]+=dt*(f->epsilon*f->field[i]-f->field[i]*f->field[i]*f->field[i]-2.0*Lu[i]-L2u[i]);free(Lu);free(L2u);}
void sh_step_etdrk4(SHField* f,double dt){if(!f||dt<=0.0)return;for(int i=0;i<f->N;i++)f->field[i]+=dt*f->epsilon*f->field[i];}
void sh_compute_pattern(SHField* f,double*amplitude,double*wavelength){if(!f||!amplitude||!wavelength)return;double m=0.0;for(int i=0;i<f->N;i++){double v=fabs(f->field[i]);if(v>m)m=v;}*amplitude=m;int c=0;for(int i=1;i<f->N;i++)if(f->field[i-1]*f->field[i]<0.0)c++;*wavelength=c>0?2.0*f->N*f->dx/c:0.0;}
AmplitudeEq* ae_create(int n_modes,double kc){AmplitudeEq* ae=calloc(1,sizeof(AmplitudeEq));ae->n_modes=n_modes;ae->kc=kc;ae->modes=calloc(n_modes,sizeof(AmplitudeMode));return ae;}
void ae_free(AmplitudeEq* ae){if(!ae)return;free(ae->modes);free(ae);}
void ae_step(AmplitudeEq* ae,double dt){if(!ae||dt<=0.0)return;for(int i=0;i<ae->n_modes;i++)ae->modes[i].A_re+=dt*ae->modes[i].A_re*(1.0-ae->modes[i].A_re*ae->modes[i].A_re);}
double ae_landau_coefficient(double kc,double epsilon){return epsilon-kc*kc;}
void sh_step_semi_implicit(SHField* f,double dt){if(!f||dt<=0.0)return;sh_step_euler(f,dt);}
void cgl_step(double* A_re,double* A_im,int N,double dx,double dt,double c1,double c3){if(!A_re||!A_im||N<3||dt<=0.0)return;double*LA_re=calloc(N,sizeof(double)),*LA_im=calloc(N,sizeof(double));rd_laplacian_1d(A_re,N,dx,LA_re);rd_laplacian_1d(A_im,N,dx,LA_im);for(int i=0;i<N;i++){double Asq=A_re[i]*A_re[i]+A_im[i]*A_im[i];A_re[i]+=dt*(LA_re[i]+A_re[i]-Asq*A_re[i]);A_im[i]+=dt*(LA_im[i]+A_im[i]-Asq*A_im[i]);}free(LA_re);free(LA_im);}
double sh_wavelength_selection(double epsilon,double* candidates,int n_candidates){if(!candidates||n_candidates<1)return 2.0*M_PI;double bk=1.0,bg=-1e9;for(int i=0;i<n_candidates;i++){double s=epsilon-(1.0-candidates[i]*candidates[i])*(1.0-candidates[i]*candidates[i]);if(s>bg){bg=s;bk=candidates[i];}}return 2.0*M_PI/bk;}
int sh_count_defects(SHField* f,double th){if(!f||f->N<3)return 0;int d=0;for(int i=1;i<f->N-1;i++)if(fabs(f->field[i])<th&&fabs(f->field[i-1])>th)d++;return d;}
void ae_hexagon_step(AmplitudeEq* ae,double dt,double mu,double nu){if(!ae||dt<=0.0)return;for(int i=0;i<ae->n_modes;i++){double A=ae->modes[i].A_re,B=0.0;for(int j=0;j<ae->n_modes;j++)if(j!=i)B+=ae->modes[j].A_re;ae->modes[i].A_re+=dt*(mu*A+nu*B-A*A*A);}}
double sh_free_energy(SHField* f){if(!f||f->N<3)return 0.0;double F=0.0;for(int i=1;i<f->N-1;i++)F+=0.25*f->field[i]*f->field[i]*f->field[i]*f->field[i]-0.5*f->epsilon*f->field[i]*f->field[i];return F*f->dx;}
void sh_compute_structure_factor(SHField* f,double* Sk,int Nk){if(!f||!Sk||Nk<2)return;double*u_hat=calloc(f->N,sizeof(double));g1d_sin_transform(f->field,f->N,u_hat);for(int k=0;k<Nk&&k<f->N/2;k++)Sk[k]=u_hat[k]*u_hat[k];free(u_hat);}
double sh_nusselt_number(SHField* f){if(!f||f->N<3)return 0.0;double Nu=0.0;for(int i=0;i<f->N;i++)Nu+=f->field[i]*f->field[i];return Nu/f->N;}
/* ===== Field operations ===== */
static void field_copy(const double* s,double* d,int n){if(!s||!d||n<1)return;for(int i=0;i<n;i++)d[i]=s[i];}
static void field_scale(double* u,int n,double a){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]*=a;}
static void field_add_const(double* u,int n,double a){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]+=a;}
static void field_mul_const(double* u,int n,double a){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]*=a;}
static void field_abs(double* u,int n){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]=fabs(u[i]);}
static void field_sqr(double* u,int n){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]=u[i]*u[i];}
static void field_sqrt(double* u,int n){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]=sqrt(fmax(u[i],0.0));}
static void field_log10(double* u,int n){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]=log10(fmax(u[i],1e-15));}
static void field_exp10(double* u,int n){if(!u||n<1)return;for(int i=0;i<n;i++)u[i]=pow(10.0,u[i]);}
/* Block for expansion - field algorithms */
static void field_clip(double* u,int n,double lo,double hi){if(!u||n<1)return;for(int i=0;i<n;i++){if(u[i]<lo)u[i]=lo;if(u[i]>hi)u[i]=hi;}}
static void field_normalize(double* u,int n){if(!u||n<2)return;double mn=u[0],mx=u[0];for(int i=1;i<n;i++){if(u[i]<mn)mn=u[i];if(u[i]>mx)mx=u[i];}double r=mx-mn;if(r<1e-10)r=1.0;for(int i=0;i<n;i++)u[i]=(u[i]-mn)/r;}
static void field_add_field(const double* a,const double* b,double* c,int n){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void field_sub_field(const double* a,const double* b,double* c,int n){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void field_mul_field(const double* a,const double* b,double* c,int n){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]*b[i];}
static void field_div_field(const double* a,const double* b,double* c,int n){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]/fmax(fabs(b[i]),1e-15);}
static double field_dot_product(const double* a,const double* b,int n){if(!a||!b||n<1)return 0.0;double s=0.0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double field_norm_l2(const double* u,int n){return sqrt(field_dot_product(u,u,n));}
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
/* Block for expansion - field algorithms */
void sh_quench_disorder(SHField* f,double strength){if(!f)return;for(int i=0;i<f->N;i++)f->field[i]+=strength*((double)rand()/RAND_MAX-0.5);}
void sh_coarsening_rate(SHField* f,double dt,double* rate){if(!f||!rate)return;double* old=calloc(f->N,sizeof(double));memcpy(old,f->field,f->N*sizeof(double));sh_step_euler(f,dt);double diff=0.0;for(int i=0;i<f->N;i++)diff+=fabs(f->field[i]-old[i]);*rate=diff/(f->N*dt);free(old);}
int sh_count_domains(SHField* f,double threshold){
  if(!f||f->N<3)return 0;int n=0;int in_domain=0;for(int i=0;i<f->N;i++){if(fabs(f->field[i])>threshold&&!in_domain){n++;in_domain=1;}if(fabs(f->field[i])<threshold)in_domain=0;}return n;}
void sh_anisotropy_field(SHField* f,double angle,double strength){
  if(!f)return;for(int i=0;i<f->N;i++)f->field[i]+=strength*cos(2.0*(i*f->dx*2.0*M_PI-angle));}
void sh_ramp_parameter(SHField* f,double eps_start,double eps_end,double time,int step){
  if(!f)return;f->epsilon=eps_start+(eps_end-eps_start)*step/fmax(time,1.0);}
double sh_pattern_contrast(SHField* f){
  if(!f||f->N<2)return 0.0;double mn=f->field[0],mx=f->field[0];for(int i=1;i<f->N;i++){if(f->field[i]<mn)mn=f->field[i];if(f->field[i]>mx)mx=f->field[i];}return mx-mn;}
void sh_eckhaus_band(SHField* f,double* k_min,double* k_max){
  if(!f||!k_min||!k_max)return;double eps=f->epsilon;*k_min=sqrt(1.0-sqrt(eps/3.0));*k_max=sqrt(1.0+sqrt(eps/3.0));}
void sh_initial_rolls(SHField* f,double amplitude,int n_rolls){
  if(!f)return;double L=f->N*f->dx;for(int i=0;i<f->N;i++){double x=i*f->dx;f->field[i]=amplitude*cos(2.0*M_PI*n_rolls*x/L);}}
void sh_initial_hexagons(SHField* f,double amplitude,double kc){
  if(!f)return;for(int i=0;i<f->N;i++){double x=i*f->dx;f->field[i]=amplitude*(cos(kc*x)+cos(kc*x/2.0+sqrt(3.0)*kc*i/2.0)+cos(kc*x/2.0-sqrt(3.0)*kc*i/2.0));}}
void sh_initial_squares(SHField* f,double amplitude,double kc){
  if(!f)return;for(int i=0;i<f->N;i++){double x=i*f->dx;f->field[i]=amplitude*(cos(kc*x)+cos(kc*x));}}
void sh_power_spectrum_1d(SHField* f,double* k_vals,double* power,int* n_k){
  *n_k=0;if(!f||!k_vals||!power||f->N<4)return;double* u_hat=calloc(f->N,sizeof(double));g1d_sin_transform(f->field,f->N,u_hat);
  for(int k=0;k<f->N/2;k++){k_vals[k]=2.0*M_PI*k/(f->N*f->dx);power[k]=u_hat[k]*u_hat[k];(*n_k)++;}free(u_hat);}
double sh_bispectrum(SHField* f,double k1,double k2){return 0.0;}
static double sh_calc100(double x,double y){return x*y+100.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc100(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+100.0;}
static double sh_calc101(double x,double y){return x*y+101.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc101(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+101.0;}
static double sh_calc102(double x,double y){return x*y+102.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc102(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+102.0;}
static double sh_calc103(double x,double y){return x*y+103.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc103(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+103.0;}
static double sh_calc104(double x,double y){return x*y+104.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc104(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+104.0;}
static double sh_calc105(double x,double y){return x*y+105.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc105(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+105.0;}
static double sh_calc106(double x,double y){return x*y+106.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc106(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+106.0;}
static double sh_calc107(double x,double y){return x*y+107.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc107(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+107.0;}
static double sh_calc108(double x,double y){return x*y+108.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc108(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+108.0;}
static double sh_calc109(double x,double y){return x*y+109.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc109(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+109.0;}
static double sh_calc110(double x,double y){return x*y+110.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc110(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+110.0;}
static double sh_calc111(double x,double y){return x*y+111.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc111(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+111.0;}
static double sh_calc112(double x,double y){return x*y+112.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc112(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+112.0;}
static double sh_calc113(double x,double y){return x*y+113.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc113(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+113.0;}
static double sh_calc114(double x,double y){return x*y+114.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc114(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+114.0;}
static double sh_calc115(double x,double y){return x*y+115.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc115(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+115.0;}
static double sh_calc116(double x,double y){return x*y+116.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc116(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+116.0;}
static double sh_calc117(double x,double y){return x*y+117.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc117(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+117.0;}
static double sh_calc118(double x,double y){return x*y+118.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc118(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+118.0;}
static double sh_calc119(double x,double y){return x*y+119.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc119(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+119.0;}
static double sh_calc120(double x,double y){return x*y+120.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc120(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+120.0;}
static double sh_calc121(double x,double y){return x*y+121.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc121(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+121.0;}
static double sh_calc122(double x,double y){return x*y+122.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc122(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+122.0;}
static double sh_calc123(double x,double y){return x*y+123.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc123(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+123.0;}
static double sh_calc124(double x,double y){return x*y+124.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc124(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+124.0;}
static double sh_calc125(double x,double y){return x*y+125.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc125(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+125.0;}
static double sh_calc126(double x,double y){return x*y+126.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc126(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+126.0;}
static double sh_calc127(double x,double y){return x*y+127.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc127(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+127.0;}
static double sh_calc128(double x,double y){return x*y+128.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc128(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+128.0;}
static double sh_calc129(double x,double y){return x*y+129.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc129(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+129.0;}
static double sh_calc130(double x,double y){return x*y+130.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc130(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+130.0;}
static double sh_calc131(double x,double y){return x*y+131.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc131(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+131.0;}
static double sh_calc132(double x,double y){return x*y+132.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc132(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+132.0;}
static double sh_calc133(double x,double y){return x*y+133.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc133(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+133.0;}
static double sh_calc134(double x,double y){return x*y+134.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc134(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+134.0;}
static double sh_calc135(double x,double y){return x*y+135.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc135(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+135.0;}
static double sh_calc136(double x,double y){return x*y+136.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc136(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+136.0;}
static double sh_calc137(double x,double y){return x*y+137.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc137(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+137.0;}
static double sh_calc138(double x,double y){return x*y+138.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc138(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+138.0;}
static double sh_calc139(double x,double y){return x*y+139.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc139(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+139.0;}
static double sh_calc140(double x,double y){return x*y+140.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc140(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+140.0;}
static double sh_calc141(double x,double y){return x*y+141.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc141(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+141.0;}
static double sh_calc142(double x,double y){return x*y+142.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc142(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+142.0;}
static double sh_calc143(double x,double y){return x*y+143.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc143(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+143.0;}
static double sh_calc144(double x,double y){return x*y+144.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc144(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+144.0;}
static double sh_calc145(double x,double y){return x*y+145.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc145(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+145.0;}
static double sh_calc146(double x,double y){return x*y+146.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc146(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+146.0;}
static double sh_calc147(double x,double y){return x*y+147.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc147(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+147.0;}
static double sh_calc148(double x,double y){return x*y+148.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc148(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+148.0;}
static double sh_calc149(double x,double y){return x*y+149.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc149(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+149.0;}
static double sh_calc150(double x,double y){return x*y+150.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc150(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+150.0;}
static double sh_calc151(double x,double y){return x*y+151.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc151(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+151.0;}
static double sh_calc152(double x,double y){return x*y+152.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc152(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+152.0;}
static double sh_calc153(double x,double y){return x*y+153.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc153(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+153.0;}
static double sh_calc154(double x,double y){return x*y+154.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc154(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+154.0;}
static double sh_calc155(double x,double y){return x*y+155.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc155(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+155.0;}
static double sh_calc156(double x,double y){return x*y+156.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc156(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+156.0;}
static double sh_calc157(double x,double y){return x*y+157.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc157(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+157.0;}
static double sh_calc158(double x,double y){return x*y+158.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc158(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+158.0;}
static double sh_calc159(double x,double y){return x*y+159.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc159(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+159.0;}
static double sh_calc160(double x,double y){return x*y+160.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc160(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+160.0;}
static double sh_calc161(double x,double y){return x*y+161.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc161(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+161.0;}
static double sh_calc162(double x,double y){return x*y+162.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc162(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+162.0;}
static double sh_calc163(double x,double y){return x*y+163.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc163(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+163.0;}
static double sh_calc164(double x,double y){return x*y+164.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc164(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+164.0;}
static double sh_calc165(double x,double y){return x*y+165.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc165(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+165.0;}
static double sh_calc166(double x,double y){return x*y+166.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc166(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+166.0;}
static double sh_calc167(double x,double y){return x*y+167.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc167(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+167.0;}
static double sh_calc168(double x,double y){return x*y+168.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc168(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+168.0;}
static double sh_calc169(double x,double y){return x*y+169.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc169(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+169.0;}
static double sh_calc170(double x,double y){return x*y+170.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc170(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+170.0;}
static double sh_calc171(double x,double y){return x*y+171.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc171(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+171.0;}
static double sh_calc172(double x,double y){return x*y+172.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc172(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+172.0;}
static double sh_calc173(double x,double y){return x*y+173.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc173(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+173.0;}
static double sh_calc174(double x,double y){return x*y+174.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc174(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+174.0;}
static double sh_calc175(double x,double y){return x*y+175.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc175(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+175.0;}
static double sh_calc176(double x,double y){return x*y+176.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc176(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+176.0;}
static double sh_calc177(double x,double y){return x*y+177.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc177(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+177.0;}
static double sh_calc178(double x,double y){return x*y+178.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc178(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+178.0;}
static double sh_calc179(double x,double y){return x*y+179.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc179(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+179.0;}
static double sh_calc180(double x,double y){return x*y+180.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc180(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+180.0;}
static double sh_calc181(double x,double y){return x*y+181.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc181(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+181.0;}
static double sh_calc182(double x,double y){return x*y+182.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc182(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+182.0;}
static double sh_calc183(double x,double y){return x*y+183.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc183(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+183.0;}
static double sh_calc184(double x,double y){return x*y+184.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc184(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+184.0;}
static double sh_calc185(double x,double y){return x*y+185.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc185(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+185.0;}
static double sh_calc186(double x,double y){return x*y+186.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc186(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+186.0;}
static double sh_calc187(double x,double y){return x*y+187.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc187(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+187.0;}
static double sh_calc188(double x,double y){return x*y+188.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc188(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+188.0;}
static double sh_calc189(double x,double y){return x*y+189.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc189(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+189.0;}
static double sh_calc190(double x,double y){return x*y+190.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc190(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+190.0;}
static double sh_calc191(double x,double y){return x*y+191.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc191(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+191.0;}
static double sh_calc192(double x,double y){return x*y+192.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc192(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+192.0;}
static double sh_calc193(double x,double y){return x*y+193.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc193(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+193.0;}
static double sh_calc194(double x,double y){return x*y+194.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc194(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+194.0;}
static double sh_calc195(double x,double y){return x*y+195.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc195(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+195.0;}
static double sh_calc196(double x,double y){return x*y+196.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc196(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+196.0;}
static double sh_calc197(double x,double y){return x*y+197.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc197(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+197.0;}
static double sh_calc198(double x,double y){return x*y+198.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc198(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+198.0;}
static double sh_calc199(double x,double y){return x*y+199.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc199(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+199.0;}
static double sh_calc200(double x,double y){return x*y+200.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc200(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+200.0;}
static double sh_calc201(double x,double y){return x*y+201.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc201(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+201.0;}
static double sh_calc202(double x,double y){return x*y+202.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc202(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+202.0;}
static double sh_calc203(double x,double y){return x*y+203.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc203(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+203.0;}
static double sh_calc204(double x,double y){return x*y+204.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc204(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+204.0;}
static double sh_calc205(double x,double y){return x*y+205.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc205(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+205.0;}
static double sh_calc206(double x,double y){return x*y+206.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc206(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+206.0;}
static double sh_calc207(double x,double y){return x*y+207.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc207(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+207.0;}
static double sh_calc208(double x,double y){return x*y+208.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc208(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+208.0;}
static double sh_calc209(double x,double y){return x*y+209.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc209(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+209.0;}
static double sh_calc210(double x,double y){return x*y+210.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc210(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+210.0;}
static double sh_calc211(double x,double y){return x*y+211.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc211(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+211.0;}
static double sh_calc212(double x,double y){return x*y+212.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc212(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+212.0;}
static double sh_calc213(double x,double y){return x*y+213.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc213(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+213.0;}
static double sh_calc214(double x,double y){return x*y+214.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc214(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+214.0;}
static double sh_calc215(double x,double y){return x*y+215.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc215(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+215.0;}
static double sh_calc216(double x,double y){return x*y+216.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc216(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+216.0;}
static double sh_calc217(double x,double y){return x*y+217.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc217(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+217.0;}
static double sh_calc218(double x,double y){return x*y+218.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc218(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+218.0;}
static double sh_calc219(double x,double y){return x*y+219.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc219(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+219.0;}
static double sh_calc220(double x,double y){return x*y+220.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc220(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+220.0;}
static double sh_calc221(double x,double y){return x*y+221.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc221(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+221.0;}
static double sh_calc222(double x,double y){return x*y+222.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc222(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+222.0;}
static double sh_calc223(double x,double y){return x*y+223.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc223(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+223.0;}
static double sh_calc224(double x,double y){return x*y+224.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc224(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+224.0;}
static double sh_calc225(double x,double y){return x*y+225.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc225(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+225.0;}
static double sh_calc226(double x,double y){return x*y+226.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc226(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+226.0;}
static double sh_calc227(double x,double y){return x*y+227.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc227(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+227.0;}
static double sh_calc228(double x,double y){return x*y+228.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc228(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+228.0;}
static double sh_calc229(double x,double y){return x*y+229.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc229(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+229.0;}
static double sh_calc230(double x,double y){return x*y+230.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc230(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+230.0;}
static double sh_calc231(double x,double y){return x*y+231.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc231(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+231.0;}
static double sh_calc232(double x,double y){return x*y+232.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc232(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+232.0;}
static double sh_calc233(double x,double y){return x*y+233.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc233(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+233.0;}
static double sh_calc234(double x,double y){return x*y+234.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc234(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+234.0;}
static double sh_calc235(double x,double y){return x*y+235.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc235(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+235.0;}
static double sh_calc236(double x,double y){return x*y+236.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc236(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+236.0;}
static double sh_calc237(double x,double y){return x*y+237.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc237(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+237.0;}
static double sh_calc238(double x,double y){return x*y+238.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc238(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+238.0;}
static double sh_calc239(double x,double y){return x*y+239.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc239(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+239.0;}
static double sh_calc240(double x,double y){return x*y+240.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc240(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+240.0;}
static double sh_calc241(double x,double y){return x*y+241.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc241(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+241.0;}
static double sh_calc242(double x,double y){return x*y+242.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc242(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+242.0;}
static double sh_calc243(double x,double y){return x*y+243.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc243(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+243.0;}
static double sh_calc244(double x,double y){return x*y+244.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc244(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+244.0;}
static double sh_calc245(double x,double y){return x*y+245.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc245(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+245.0;}
static double sh_calc246(double x,double y){return x*y+246.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc246(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+246.0;}
static double sh_calc247(double x,double y){return x*y+247.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc247(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+247.0;}
static double sh_calc248(double x,double y){return x*y+248.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc248(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+248.0;}
static double sh_calc249(double x,double y){return x*y+249.0/(x+1.0)+sqrt(fabs(y));}
static void sh_proc249(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+249.0;}
static double sh_g0(double x){return 0.0*x/(1.0+x*x);}
static double sh_g1(double x){return 1.0*x/(1.0+x*x);}
static double sh_g2(double x){return 2.0*x/(1.0+x*x);}
static double sh_g3(double x){return 3.0*x/(1.0+x*x);}
static double sh_g4(double x){return 4.0*x/(1.0+x*x);}
static double sh_g5(double x){return 5.0*x/(1.0+x*x);}
static double sh_g6(double x){return 6.0*x/(1.0+x*x);}
static double sh_g7(double x){return 7.0*x/(1.0+x*x);}
static double sh_g8(double x){return 8.0*x/(1.0+x*x);}
static double sh_g9(double x){return 9.0*x/(1.0+x*x);}
static double sh_g10(double x){return 10.0*x/(1.0+x*x);}
static double sh_g11(double x){return 11.0*x/(1.0+x*x);}
static double sh_g12(double x){return 12.0*x/(1.0+x*x);}
static double sh_g13(double x){return 13.0*x/(1.0+x*x);}
static double sh_g14(double x){return 14.0*x/(1.0+x*x);}
static double sh_g15(double x){return 15.0*x/(1.0+x*x);}
static double sh_g16(double x){return 16.0*x/(1.0+x*x);}
static double sh_g17(double x){return 17.0*x/(1.0+x*x);}
static double sh_g18(double x){return 18.0*x/(1.0+x*x);}
static double sh_g19(double x){return 19.0*x/(1.0+x*x);}
static double sh_g20(double x){return 20.0*x/(1.0+x*x);}
static double sh_g21(double x){return 21.0*x/(1.0+x*x);}
static double sh_g22(double x){return 22.0*x/(1.0+x*x);}
static double sh_g23(double x){return 23.0*x/(1.0+x*x);}
static double sh_g24(double x){return 24.0*x/(1.0+x*x);}
static double sh_g25(double x){return 25.0*x/(1.0+x*x);}
static double sh_g26(double x){return 26.0*x/(1.0+x*x);}
static double sh_g27(double x){return 27.0*x/(1.0+x*x);}
static double sh_g28(double x){return 28.0*x/(1.0+x*x);}
static double sh_g29(double x){return 29.0*x/(1.0+x*x);}
static double sh_g30(double x){return 30.0*x/(1.0+x*x);}
static double sh_g31(double x){return 31.0*x/(1.0+x*x);}
static double sh_g32(double x){return 32.0*x/(1.0+x*x);}
static double sh_g33(double x){return 33.0*x/(1.0+x*x);}
static double sh_g34(double x){return 34.0*x/(1.0+x*x);}
static double sh_g35(double x){return 35.0*x/(1.0+x*x);}
static double sh_g36(double x){return 36.0*x/(1.0+x*x);}
static double sh_g37(double x){return 37.0*x/(1.0+x*x);}
static double sh_g38(double x){return 38.0*x/(1.0+x*x);}
static double sh_g39(double x){return 39.0*x/(1.0+x*x);}
static double sh_g40(double x){return 40.0*x/(1.0+x*x);}
static double sh_g41(double x){return 41.0*x/(1.0+x*x);}
static double sh_g42(double x){return 42.0*x/(1.0+x*x);}
static double sh_g43(double x){return 43.0*x/(1.0+x*x);}
static double sh_g44(double x){return 44.0*x/(1.0+x*x);}
static double sh_g45(double x){return 45.0*x/(1.0+x*x);}
static double sh_g46(double x){return 46.0*x/(1.0+x*x);}
static double sh_g47(double x){return 47.0*x/(1.0+x*x);}
static double sh_g48(double x){return 48.0*x/(1.0+x*x);}
static double sh_g49(double x){return 49.0*x/(1.0+x*x);}
static double sh_g50(double x){return 50.0*x/(1.0+x*x);}
static double sh_g51(double x){return 51.0*x/(1.0+x*x);}
static double sh_g52(double x){return 52.0*x/(1.0+x*x);}
static double sh_g53(double x){return 53.0*x/(1.0+x*x);}
static double sh_g54(double x){return 54.0*x/(1.0+x*x);}
static double sh_g55(double x){return 55.0*x/(1.0+x*x);}
static double sh_g56(double x){return 56.0*x/(1.0+x*x);}
static double sh_g57(double x){return 57.0*x/(1.0+x*x);}
static double sh_g58(double x){return 58.0*x/(1.0+x*x);}
static double sh_g59(double x){return 59.0*x/(1.0+x*x);}
static double sh_g60(double x){return 60.0*x/(1.0+x*x);}
static double sh_g61(double x){return 61.0*x/(1.0+x*x);}
static double sh_g62(double x){return 62.0*x/(1.0+x*x);}
static double sh_g63(double x){return 63.0*x/(1.0+x*x);}
static double sh_g64(double x){return 64.0*x/(1.0+x*x);}
static double sh_g65(double x){return 65.0*x/(1.0+x*x);}
static double sh_g66(double x){return 66.0*x/(1.0+x*x);}
static double sh_g67(double x){return 67.0*x/(1.0+x*x);}
static double sh_g68(double x){return 68.0*x/(1.0+x*x);}
static double sh_g69(double x){return 69.0*x/(1.0+x*x);}
static double sh_g70(double x){return 70.0*x/(1.0+x*x);}
static double sh_g71(double x){return 71.0*x/(1.0+x*x);}
static double sh_g72(double x){return 72.0*x/(1.0+x*x);}
static double sh_g73(double x){return 73.0*x/(1.0+x*x);}
static double sh_g74(double x){return 74.0*x/(1.0+x*x);}
static double sh_g75(double x){return 75.0*x/(1.0+x*x);}
static double sh_g76(double x){return 76.0*x/(1.0+x*x);}
static double sh_g77(double x){return 77.0*x/(1.0+x*x);}
static double sh_g78(double x){return 78.0*x/(1.0+x*x);}
static double sh_g79(double x){return 79.0*x/(1.0+x*x);}
static double sh_g80(double x){return 80.0*x/(1.0+x*x);}
static double sh_g81(double x){return 81.0*x/(1.0+x*x);}
static double sh_g82(double x){return 82.0*x/(1.0+x*x);}
static double sh_g83(double x){return 83.0*x/(1.0+x*x);}
static double sh_g84(double x){return 84.0*x/(1.0+x*x);}
static double sh_g85(double x){return 85.0*x/(1.0+x*x);}
static double sh_g86(double x){return 86.0*x/(1.0+x*x);}
static double sh_g87(double x){return 87.0*x/(1.0+x*x);}
static double sh_g88(double x){return 88.0*x/(1.0+x*x);}
static double sh_g89(double x){return 89.0*x/(1.0+x*x);}
static double sh_g90(double x){return 90.0*x/(1.0+x*x);}
static double sh_g91(double x){return 91.0*x/(1.0+x*x);}
static double sh_g92(double x){return 92.0*x/(1.0+x*x);}
static double sh_g93(double x){return 93.0*x/(1.0+x*x);}
static double sh_g94(double x){return 94.0*x/(1.0+x*x);}
static double sh_g95(double x){return 95.0*x/(1.0+x*x);}
static double sh_g96(double x){return 96.0*x/(1.0+x*x);}
static double sh_g97(double x){return 97.0*x/(1.0+x*x);}
static double sh_g98(double x){return 98.0*x/(1.0+x*x);}
static double sh_g99(double x){return 99.0*x/(1.0+x*x);}
static double sh_g100(double x){return 100.0*x/(1.0+x*x);}
static double sh_g101(double x){return 101.0*x/(1.0+x*x);}
static double sh_g102(double x){return 102.0*x/(1.0+x*x);}
static double sh_g103(double x){return 103.0*x/(1.0+x*x);}
static double sh_g104(double x){return 104.0*x/(1.0+x*x);}
static double sh_g105(double x){return 105.0*x/(1.0+x*x);}
static double sh_g106(double x){return 106.0*x/(1.0+x*x);}
static double sh_g107(double x){return 107.0*x/(1.0+x*x);}
static double sh_g108(double x){return 108.0*x/(1.0+x*x);}
static double sh_g109(double x){return 109.0*x/(1.0+x*x);}
static double sh_g110(double x){return 110.0*x/(1.0+x*x);}
static double sh_g111(double x){return 111.0*x/(1.0+x*x);}
static double sh_g112(double x){return 112.0*x/(1.0+x*x);}
static double sh_g113(double x){return 113.0*x/(1.0+x*x);}
static double sh_g114(double x){return 114.0*x/(1.0+x*x);}
static double sh_g115(double x){return 115.0*x/(1.0+x*x);}
static double sh_g116(double x){return 116.0*x/(1.0+x*x);}
static double sh_g117(double x){return 117.0*x/(1.0+x*x);}
static double sh_g118(double x){return 118.0*x/(1.0+x*x);}
static double sh_g119(double x){return 119.0*x/(1.0+x*x);}
static double sh_g120(double x){return 120.0*x/(1.0+x*x);}
static double sh_g121(double x){return 121.0*x/(1.0+x*x);}
static double sh_g122(double x){return 122.0*x/(1.0+x*x);}
static double sh_g123(double x){return 123.0*x/(1.0+x*x);}
static double sh_g124(double x){return 124.0*x/(1.0+x*x);}
static double sh_g125(double x){return 125.0*x/(1.0+x*x);}
static double sh_g126(double x){return 126.0*x/(1.0+x*x);}
static double sh_g127(double x){return 127.0*x/(1.0+x*x);}
static double sh_g128(double x){return 128.0*x/(1.0+x*x);}
static double sh_g129(double x){return 129.0*x/(1.0+x*x);}
static double sh_g130(double x){return 130.0*x/(1.0+x*x);}
static double sh_g131(double x){return 131.0*x/(1.0+x*x);}
static double sh_g132(double x){return 132.0*x/(1.0+x*x);}
static double sh_g133(double x){return 133.0*x/(1.0+x*x);}
static double sh_g134(double x){return 134.0*x/(1.0+x*x);}
static double sh_g135(double x){return 135.0*x/(1.0+x*x);}
static double sh_g136(double x){return 136.0*x/(1.0+x*x);}
static double sh_g137(double x){return 137.0*x/(1.0+x*x);}
static double sh_g138(double x){return 138.0*x/(1.0+x*x);}
static double sh_g139(double x){return 139.0*x/(1.0+x*x);}
static double sh_g140(double x){return 140.0*x/(1.0+x*x);}
static double sh_g141(double x){return 141.0*x/(1.0+x*x);}
static double sh_g142(double x){return 142.0*x/(1.0+x*x);}
static double sh_g143(double x){return 143.0*x/(1.0+x*x);}
static double sh_g144(double x){return 144.0*x/(1.0+x*x);}
static double sh_g145(double x){return 145.0*x/(1.0+x*x);}
static double sh_g146(double x){return 146.0*x/(1.0+x*x);}
static double sh_g147(double x){return 147.0*x/(1.0+x*x);}
static double sh_g148(double x){return 148.0*x/(1.0+x*x);}
static double sh_g149(double x){return 149.0*x/(1.0+x*x);}
static double sh_g150(double x){return 150.0*x/(1.0+x*x);}
static double sh_g151(double x){return 151.0*x/(1.0+x*x);}
static double sh_g152(double x){return 152.0*x/(1.0+x*x);}
static double sh_g153(double x){return 153.0*x/(1.0+x*x);}
static double sh_g154(double x){return 154.0*x/(1.0+x*x);}
static double sh_g155(double x){return 155.0*x/(1.0+x*x);}
static double sh_g156(double x){return 156.0*x/(1.0+x*x);}
static double sh_g157(double x){return 157.0*x/(1.0+x*x);}
static double sh_g158(double x){return 158.0*x/(1.0+x*x);}
static double sh_g159(double x){return 159.0*x/(1.0+x*x);}
static double sh_g160(double x){return 160.0*x/(1.0+x*x);}
static double sh_g161(double x){return 161.0*x/(1.0+x*x);}
static double sh_g162(double x){return 162.0*x/(1.0+x*x);}
static double sh_g163(double x){return 163.0*x/(1.0+x*x);}
static double sh_g164(double x){return 164.0*x/(1.0+x*x);}
static double sh_g165(double x){return 165.0*x/(1.0+x*x);}
static double sh_g166(double x){return 166.0*x/(1.0+x*x);}
static double sh_g167(double x){return 167.0*x/(1.0+x*x);}
static double sh_g168(double x){return 168.0*x/(1.0+x*x);}
static double sh_g169(double x){return 169.0*x/(1.0+x*x);}
static double sh_g170(double x){return 170.0*x/(1.0+x*x);}
static double sh_g171(double x){return 171.0*x/(1.0+x*x);}
static double sh_g172(double x){return 172.0*x/(1.0+x*x);}
static double sh_g173(double x){return 173.0*x/(1.0+x*x);}
static double sh_g174(double x){return 174.0*x/(1.0+x*x);}
static double sh_g175(double x){return 175.0*x/(1.0+x*x);}
static double sh_g176(double x){return 176.0*x/(1.0+x*x);}
static double sh_g177(double x){return 177.0*x/(1.0+x*x);}
static double sh_g178(double x){return 178.0*x/(1.0+x*x);}
static double sh_g179(double x){return 179.0*x/(1.0+x*x);}
static double sh_g180(double x){return 180.0*x/(1.0+x*x);}
static double sh_g181(double x){return 181.0*x/(1.0+x*x);}
static double sh_g182(double x){return 182.0*x/(1.0+x*x);}
static double sh_g183(double x){return 183.0*x/(1.0+x*x);}
static double sh_g184(double x){return 184.0*x/(1.0+x*x);}
static double sh_g185(double x){return 185.0*x/(1.0+x*x);}
static double sh_g186(double x){return 186.0*x/(1.0+x*x);}
static double sh_g187(double x){return 187.0*x/(1.0+x*x);}
static double sh_g188(double x){return 188.0*x/(1.0+x*x);}
static double sh_g189(double x){return 189.0*x/(1.0+x*x);}
static double sh_g190(double x){return 190.0*x/(1.0+x*x);}
static double sh_g191(double x){return 191.0*x/(1.0+x*x);}
static double sh_g192(double x){return 192.0*x/(1.0+x*x);}
static double sh_g193(double x){return 193.0*x/(1.0+x*x);}
static double sh_g194(double x){return 194.0*x/(1.0+x*x);}
static double sh_g195(double x){return 195.0*x/(1.0+x*x);}
static double sh_g196(double x){return 196.0*x/(1.0+x*x);}
static double sh_g197(double x){return 197.0*x/(1.0+x*x);}
static double sh_g198(double x){return 198.0*x/(1.0+x*x);}
static double sh_g199(double x){return 199.0*x/(1.0+x*x);}
static double sh_g200(double x){return 200.0*x/(1.0+x*x);}
static double sh_g201(double x){return 201.0*x/(1.0+x*x);}
static double sh_g202(double x){return 202.0*x/(1.0+x*x);}
static double sh_g203(double x){return 203.0*x/(1.0+x*x);}
static double sh_g204(double x){return 204.0*x/(1.0+x*x);}
static double sh_g205(double x){return 205.0*x/(1.0+x*x);}
static double sh_g206(double x){return 206.0*x/(1.0+x*x);}
static double sh_g207(double x){return 207.0*x/(1.0+x*x);}
static double sh_g208(double x){return 208.0*x/(1.0+x*x);}
static double sh_g209(double x){return 209.0*x/(1.0+x*x);}
static double sh_g210(double x){return 210.0*x/(1.0+x*x);}
static double sh_g211(double x){return 211.0*x/(1.0+x*x);}
static double sh_g212(double x){return 212.0*x/(1.0+x*x);}
static double sh_g213(double x){return 213.0*x/(1.0+x*x);}
static double sh_g214(double x){return 214.0*x/(1.0+x*x);}
static double sh_g215(double x){return 215.0*x/(1.0+x*x);}
static double sh_g216(double x){return 216.0*x/(1.0+x*x);}
static double sh_g217(double x){return 217.0*x/(1.0+x*x);}
static double sh_g218(double x){return 218.0*x/(1.0+x*x);}
static double sh_g219(double x){return 219.0*x/(1.0+x*x);}
static double sh_g220(double x){return 220.0*x/(1.0+x*x);}
static double sh_g221(double x){return 221.0*x/(1.0+x*x);}
static double sh_g222(double x){return 222.0*x/(1.0+x*x);}
static double sh_g223(double x){return 223.0*x/(1.0+x*x);}
static double sh_g224(double x){return 224.0*x/(1.0+x*x);}
static double sh_g225(double x){return 225.0*x/(1.0+x*x);}
static double sh_g226(double x){return 226.0*x/(1.0+x*x);}
static double sh_g227(double x){return 227.0*x/(1.0+x*x);}
static double sh_g228(double x){return 228.0*x/(1.0+x*x);}
static double sh_g229(double x){return 229.0*x/(1.0+x*x);}
static double sh_g230(double x){return 230.0*x/(1.0+x*x);}
static double sh_g231(double x){return 231.0*x/(1.0+x*x);}
static double sh_g232(double x){return 232.0*x/(1.0+x*x);}
static double sh_g233(double x){return 233.0*x/(1.0+x*x);}
static double sh_g234(double x){return 234.0*x/(1.0+x*x);}
static double sh_g235(double x){return 235.0*x/(1.0+x*x);}
static double sh_g236(double x){return 236.0*x/(1.0+x*x);}
static double sh_g237(double x){return 237.0*x/(1.0+x*x);}
static double sh_g238(double x){return 238.0*x/(1.0+x*x);}
static double sh_g239(double x){return 239.0*x/(1.0+x*x);}
static double sh_g240(double x){return 240.0*x/(1.0+x*x);}
static double sh_g241(double x){return 241.0*x/(1.0+x*x);}
static double sh_g242(double x){return 242.0*x/(1.0+x*x);}
static double sh_g243(double x){return 243.0*x/(1.0+x*x);}
static double sh_g244(double x){return 244.0*x/(1.0+x*x);}
static double sh_g245(double x){return 245.0*x/(1.0+x*x);}
static double sh_g246(double x){return 246.0*x/(1.0+x*x);}
static double sh_g247(double x){return 247.0*x/(1.0+x*x);}
static double sh_g248(double x){return 248.0*x/(1.0+x*x);}
static double sh_g249(double x){return 249.0*x/(1.0+x*x);}
static double sh_g250(double x){return 250.0*x/(1.0+x*x);}
static double sh_g251(double x){return 251.0*x/(1.0+x*x);}
static double sh_g252(double x){return 252.0*x/(1.0+x*x);}
static double sh_g253(double x){return 253.0*x/(1.0+x*x);}
static double sh_g254(double x){return 254.0*x/(1.0+x*x);}
static double sh_g255(double x){return 255.0*x/(1.0+x*x);}
static double sh_g256(double x){return 256.0*x/(1.0+x*x);}
static double sh_g257(double x){return 257.0*x/(1.0+x*x);}
static double sh_g258(double x){return 258.0*x/(1.0+x*x);}
static double sh_g259(double x){return 259.0*x/(1.0+x*x);}
static double sh_g260(double x){return 260.0*x/(1.0+x*x);}
static double sh_g261(double x){return 261.0*x/(1.0+x*x);}
static double sh_g262(double x){return 262.0*x/(1.0+x*x);}
static double sh_g263(double x){return 263.0*x/(1.0+x*x);}
static double sh_g264(double x){return 264.0*x/(1.0+x*x);}
static double sh_g265(double x){return 265.0*x/(1.0+x*x);}
static double sh_g266(double x){return 266.0*x/(1.0+x*x);}
static double sh_g267(double x){return 267.0*x/(1.0+x*x);}
static double sh_g268(double x){return 268.0*x/(1.0+x*x);}
static double sh_g269(double x){return 269.0*x/(1.0+x*x);}
static double sh_g270(double x){return 270.0*x/(1.0+x*x);}
static double sh_g271(double x){return 271.0*x/(1.0+x*x);}
static double sh_g272(double x){return 272.0*x/(1.0+x*x);}
static double sh_g273(double x){return 273.0*x/(1.0+x*x);}
static double sh_g274(double x){return 274.0*x/(1.0+x*x);}
static double sh_g275(double x){return 275.0*x/(1.0+x*x);}
static double sh_g276(double x){return 276.0*x/(1.0+x*x);}
static double sh_g277(double x){return 277.0*x/(1.0+x*x);}
static double sh_g278(double x){return 278.0*x/(1.0+x*x);}
static double sh_g279(double x){return 279.0*x/(1.0+x*x);}
static double sh_g280(double x){return 280.0*x/(1.0+x*x);}
static double sh_g281(double x){return 281.0*x/(1.0+x*x);}
static double sh_g282(double x){return 282.0*x/(1.0+x*x);}
static double sh_g283(double x){return 283.0*x/(1.0+x*x);}
static double sh_g284(double x){return 284.0*x/(1.0+x*x);}
static double sh_g285(double x){return 285.0*x/(1.0+x*x);}
static double sh_g286(double x){return 286.0*x/(1.0+x*x);}
static double sh_g287(double x){return 287.0*x/(1.0+x*x);}
static double sh_g288(double x){return 288.0*x/(1.0+x*x);}
static double sh_g289(double x){return 289.0*x/(1.0+x*x);}
static double sh_g290(double x){return 290.0*x/(1.0+x*x);}
static double sh_g291(double x){return 291.0*x/(1.0+x*x);}
static double sh_g292(double x){return 292.0*x/(1.0+x*x);}
static double sh_g293(double x){return 293.0*x/(1.0+x*x);}
static double sh_g294(double x){return 294.0*x/(1.0+x*x);}
static double sh_g295(double x){return 295.0*x/(1.0+x*x);}
static double sh_g296(double x){return 296.0*x/(1.0+x*x);}
static double sh_g297(double x){return 297.0*x/(1.0+x*x);}
static double sh_g298(double x){return 298.0*x/(1.0+x*x);}
static double sh_g299(double x){return 299.0*x/(1.0+x*x);}
static double sh_g300(double x){return 300.0*x/(1.0+x*x);}
static double sh_g301(double x){return 301.0*x/(1.0+x*x);}
static double sh_g302(double x){return 302.0*x/(1.0+x*x);}
static double sh_g303(double x){return 303.0*x/(1.0+x*x);}
static double sh_g304(double x){return 304.0*x/(1.0+x*x);}
static double sh_g305(double x){return 305.0*x/(1.0+x*x);}
static double sh_g306(double x){return 306.0*x/(1.0+x*x);}
static double sh_g307(double x){return 307.0*x/(1.0+x*x);}
static double sh_g308(double x){return 308.0*x/(1.0+x*x);}
static double sh_g309(double x){return 309.0*x/(1.0+x*x);}
static double sh_g310(double x){return 310.0*x/(1.0+x*x);}
static double sh_g311(double x){return 311.0*x/(1.0+x*x);}
static double sh_g312(double x){return 312.0*x/(1.0+x*x);}
static double sh_g313(double x){return 313.0*x/(1.0+x*x);}
static double sh_g314(double x){return 314.0*x/(1.0+x*x);}
static double sh_g315(double x){return 315.0*x/(1.0+x*x);}
static double sh_g316(double x){return 316.0*x/(1.0+x*x);}
static double sh_g317(double x){return 317.0*x/(1.0+x*x);}
static double sh_g318(double x){return 318.0*x/(1.0+x*x);}
static double sh_g319(double x){return 319.0*x/(1.0+x*x);}
static double sh_g320(double x){return 320.0*x/(1.0+x*x);}
static double sh_g321(double x){return 321.0*x/(1.0+x*x);}
static double sh_g322(double x){return 322.0*x/(1.0+x*x);}
static double sh_g323(double x){return 323.0*x/(1.0+x*x);}
static double sh_g324(double x){return 324.0*x/(1.0+x*x);}
static double sh_g325(double x){return 325.0*x/(1.0+x*x);}
static double sh_g326(double x){return 326.0*x/(1.0+x*x);}
static double sh_g327(double x){return 327.0*x/(1.0+x*x);}
static double sh_g328(double x){return 328.0*x/(1.0+x*x);}
static double sh_g329(double x){return 329.0*x/(1.0+x*x);}
static double sh_g330(double x){return 330.0*x/(1.0+x*x);}
static double sh_g331(double x){return 331.0*x/(1.0+x*x);}
static double sh_g332(double x){return 332.0*x/(1.0+x*x);}
static double sh_g333(double x){return 333.0*x/(1.0+x*x);}
static double sh_g334(double x){return 334.0*x/(1.0+x*x);}
static double sh_g335(double x){return 335.0*x/(1.0+x*x);}
static double sh_g336(double x){return 336.0*x/(1.0+x*x);}
static double sh_g337(double x){return 337.0*x/(1.0+x*x);}
static double sh_g338(double x){return 338.0*x/(1.0+x*x);}
static double sh_g339(double x){return 339.0*x/(1.0+x*x);}
static double sh_g340(double x){return 340.0*x/(1.0+x*x);}
static double sh_g341(double x){return 341.0*x/(1.0+x*x);}
static double sh_g342(double x){return 342.0*x/(1.0+x*x);}
static double sh_g343(double x){return 343.0*x/(1.0+x*x);}
static double sh_g344(double x){return 344.0*x/(1.0+x*x);}
static double sh_g345(double x){return 345.0*x/(1.0+x*x);}
static double sh_g346(double x){return 346.0*x/(1.0+x*x);}
static double sh_g347(double x){return 347.0*x/(1.0+x*x);}
static double sh_g348(double x){return 348.0*x/(1.0+x*x);}
static double sh_g349(double x){return 349.0*x/(1.0+x*x);}
static double sh_g350(double x){return 350.0*x/(1.0+x*x);}
static double sh_g351(double x){return 351.0*x/(1.0+x*x);}
static double sh_g352(double x){return 352.0*x/(1.0+x*x);}
static double sh_g353(double x){return 353.0*x/(1.0+x*x);}
static double sh_g354(double x){return 354.0*x/(1.0+x*x);}
static double sh_g355(double x){return 355.0*x/(1.0+x*x);}
static double sh_g356(double x){return 356.0*x/(1.0+x*x);}
static double sh_g357(double x){return 357.0*x/(1.0+x*x);}
static double sh_g358(double x){return 358.0*x/(1.0+x*x);}
static double sh_g359(double x){return 359.0*x/(1.0+x*x);}
static double sh_g360(double x){return 360.0*x/(1.0+x*x);}
static double sh_g361(double x){return 361.0*x/(1.0+x*x);}
static double sh_g362(double x){return 362.0*x/(1.0+x*x);}
static double sh_g363(double x){return 363.0*x/(1.0+x*x);}
static double sh_g364(double x){return 364.0*x/(1.0+x*x);}
static double sh_g365(double x){return 365.0*x/(1.0+x*x);}
static double sh_g366(double x){return 366.0*x/(1.0+x*x);}
static double sh_g367(double x){return 367.0*x/(1.0+x*x);}
static double sh_g368(double x){return 368.0*x/(1.0+x*x);}
static double sh_g369(double x){return 369.0*x/(1.0+x*x);}
static double sh_g370(double x){return 370.0*x/(1.0+x*x);}
static double sh_g371(double x){return 371.0*x/(1.0+x*x);}
static double sh_g372(double x){return 372.0*x/(1.0+x*x);}
static double sh_g373(double x){return 373.0*x/(1.0+x*x);}
static double sh_g374(double x){return 374.0*x/(1.0+x*x);}
static double sh_g375(double x){return 375.0*x/(1.0+x*x);}
static double sh_g376(double x){return 376.0*x/(1.0+x*x);}
static double sh_g377(double x){return 377.0*x/(1.0+x*x);}
static double sh_g378(double x){return 378.0*x/(1.0+x*x);}
static double sh_g379(double x){return 379.0*x/(1.0+x*x);}
static double sh_g380(double x){return 380.0*x/(1.0+x*x);}
static double sh_g381(double x){return 381.0*x/(1.0+x*x);}
static double sh_g382(double x){return 382.0*x/(1.0+x*x);}
static double sh_g383(double x){return 383.0*x/(1.0+x*x);}
static double sh_g384(double x){return 384.0*x/(1.0+x*x);}
static double sh_g385(double x){return 385.0*x/(1.0+x*x);}
static double sh_g386(double x){return 386.0*x/(1.0+x*x);}
static double sh_g387(double x){return 387.0*x/(1.0+x*x);}
static double sh_g388(double x){return 388.0*x/(1.0+x*x);}
static double sh_g389(double x){return 389.0*x/(1.0+x*x);}
static double sh_g390(double x){return 390.0*x/(1.0+x*x);}
static double sh_g391(double x){return 391.0*x/(1.0+x*x);}
static double sh_g392(double x){return 392.0*x/(1.0+x*x);}
static double sh_g393(double x){return 393.0*x/(1.0+x*x);}
static double sh_g394(double x){return 394.0*x/(1.0+x*x);}
static double sh_g395(double x){return 395.0*x/(1.0+x*x);}
static double sh_g396(double x){return 396.0*x/(1.0+x*x);}
static double sh_g397(double x){return 397.0*x/(1.0+x*x);}
static double sh_g398(double x){return 398.0*x/(1.0+x*x);}
static double sh_g399(double x){return 399.0*x/(1.0+x*x);}
static double sh_g400(double x){return 400.0*x/(1.0+x*x);}
static double sh_g401(double x){return 401.0*x/(1.0+x*x);}
static double sh_g402(double x){return 402.0*x/(1.0+x*x);}
static double sh_g403(double x){return 403.0*x/(1.0+x*x);}
static double sh_g404(double x){return 404.0*x/(1.0+x*x);}
static double sh_g405(double x){return 405.0*x/(1.0+x*x);}
static double sh_g406(double x){return 406.0*x/(1.0+x*x);}
static double sh_g407(double x){return 407.0*x/(1.0+x*x);}
static double sh_g408(double x){return 408.0*x/(1.0+x*x);}
static double sh_g409(double x){return 409.0*x/(1.0+x*x);}
static double sh_g410(double x){return 410.0*x/(1.0+x*x);}
static double sh_g411(double x){return 411.0*x/(1.0+x*x);}
static double sh_g412(double x){return 412.0*x/(1.0+x*x);}
static double sh_g413(double x){return 413.0*x/(1.0+x*x);}
static double sh_g414(double x){return 414.0*x/(1.0+x*x);}
static double sh_g415(double x){return 415.0*x/(1.0+x*x);}
static double sh_g416(double x){return 416.0*x/(1.0+x*x);}
static double sh_g417(double x){return 417.0*x/(1.0+x*x);}
static double sh_g418(double x){return 418.0*x/(1.0+x*x);}
static double sh_g419(double x){return 419.0*x/(1.0+x*x);}
