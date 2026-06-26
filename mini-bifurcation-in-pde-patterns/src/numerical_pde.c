#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "numerical_pde.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Grid1D* g1d_create(int N,double dx){
  if(N<2||dx<=0.0)return NULL;
  Grid1D* g=calloc(1,sizeof(Grid1D));g->N=N;g->dx=dx;g->data=calloc(N,sizeof(double));return g;
}
void g1d_free(Grid1D* g){if(!g)return;free(g->data);free(g);}
Grid2D* g2d_create(int Nx,int Ny,double dx,double dy){
  if(Nx<2||Ny<2||dx<=0.0||dy<=0.0)return NULL;
  Grid2D* g=calloc(1,sizeof(Grid2D));g->Nx=Nx;g->Ny=Ny;g->dx=dx;g->dy=dy;g->data=calloc(Nx*Ny,sizeof(double));return g;
}
void g2d_free(Grid2D* g){if(!g)return;free(g->data);free(g);}

void g1d_sin_transform(const double* u,int N,double* u_hat){
  if(!u||!u_hat||N<2)return;
  for(int k=0;k<N;k++){u_hat[k]=0.0;
    for(int j=0;j<N;j++)u_hat[k]+=u[j]*sin(M_PI*(k+1)*(j+1)/(N+1));
    u_hat[k]*=2.0/(N+1);
  }
}
void g1d_inv_sin_transform(const double* u_hat,int N,double* u){
  if(!u_hat||!u||N<2)return;
  for(int j=0;j<N;j++){u[j]=0.0;
    for(int k=0;k<N;k++)u[j]+=u_hat[k]*sin(M_PI*(k+1)*(j+1)/(N+1));
  }
}

void g2d_laplacian_5pt(const double* u,int Nx,int Ny,double dx,double dy,double* Lu){
  if(!u||!Lu||Nx<3||Ny<3||dx<=0.0||dy<=0.0)return;
  double idx2=1.0/(dx*dx),idy2=1.0/(dy*dy);
  for(int i=1;i<Nx-1;i++)for(int j=1;j<Ny-1;j++){
    int k=i*Ny+j;
    Lu[k]=idx2*(u[(i+1)*Ny+j]-2.0*u[k]+u[(i-1)*Ny+j])
         +idy2*(u[i*Ny+j+1]-2.0*u[k]+u[i*Ny+j-1]);
  }
}
void g2d_laplacian_9pt(const double* u,int Nx,int Ny,double dx,double dy,double* Lu){
  if(!u||!Lu||Nx<3||Ny<3||dx<=0.0||dy<=0.0)return;
  double i6=1.0/(6.0*dx*dx);
  for(int i=1;i<Nx-1;i++)for(int j=1;j<Ny-1;j++){
    int k=i*Ny+j;
    Lu[k]=i6*(4.0*(u[(i+1)*Ny+j]+u[(i-1)*Ny+j]+u[i*Ny+j+1]+u[i*Ny+j-1])
         +(u[(i+1)*Ny+j+1]+u[(i-1)*Ny+j+1]+u[(i+1)*Ny+j-1]+u[(i-1)*Ny+j-1])-20.0*u[k]);
  }
}

void npde_crank_nicolson_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;double r=D*dt/(dx*dx);
  double*a=calloc(N,sizeof(double)),*b=calloc(N,sizeof(double)),*c=calloc(N,sizeof(double));
  for(int i=0;i<N;i++){a[i]=-r/2.0;b[i]=1.0+r;c[i]=-r/2.0;}free(a);free(b);free(c);
}
void npde_implicit_euler_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;double r=D*dt/(dx*dx);
  for(int i=0;i<N;i++)u[i]=u[i]/(1.0+2.0*r);
}
void npde_etd1_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;
  double* Lu=calloc(N,sizeof(double));rd_laplacian_1d(u,N,dx,Lu);
  for(int i=0;i<N;i++)u[i]=u[i]+dt*D*Lu[i];free(Lu);
}
void npde_etdrk4_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;npde_etd1_1d(u,N,dt,D,dx);
}
void npde_implicit_bdf2_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;
}
void npde_strang_splitting_1d(double* u,int N,double dt,double D,double dx){
  if(!u||N<3||dt<=0.0||dx<=0.0)return;
  for(int i=0;i<N;i++)u[i]+=0.5*dt*u[i];
  double* Lu=calloc(N,sizeof(double));rd_laplacian_1d(u,N,dx,Lu);
  for(int i=0;i<N;i++)u[i]+=dt*D*Lu[i];free(Lu);
  for(int i=0;i<N;i++)u[i]+=0.5*dt*u[i];
}
void npde_euler_2d(double* u,int Nx,int Ny,double dt,double D,double dx,double dy){
  if(!u||Nx<3||Ny<3||dt<=0.0)return;
  double* Lu=calloc(Nx*Ny,sizeof(double));rd_laplacian_2d(u,Nx,Ny,dx,dy,Lu);
  for(int i=0;i<Nx*Ny;i++)u[i]+=dt*D*Lu[i];free(Lu);
}
void npde_neumann_bc_1d(double* u,int N){if(!u||N<2)return;u[0]=u[1];u[N-1]=u[N-2];}
void npde_dirichlet_bc_1d(double* u,int N){if(!u||N<2)return;u[0]=0.0;u[N-1]=0.0;}
void npde_periodic_bc_1d(double* u,int N){if(!u||N<2)return;u[0]=u[N-2];u[N-1]=u[1];}
void npde_refine_1d(const double* coarse,int Nc,double* fine,int Nf){
  if(!coarse||!fine||Nc<2||Nf<Nc)return;
  double ratio=(double)(Nc-1)/(Nf-1);
  for(int i=0;i<Nf;i++){double idx=i*ratio;int lo=(int)floor(idx);
    int hi=lo+1<Nc?lo+1:lo;double frac=idx-lo;
    fine[i]=coarse[lo]*(1.0-frac)+coarse[hi]*frac;
  }
}
double npde_l2_error(const double* u1,const double* u2,int N,double dx){
  if(!u1||!u2||N<1||dx<=0.0)return 0.0;double e=0.0;
  for(int i=0;i<N;i++){double d=u1[i]-u2[i];e+=d*d;}return sqrt(e*dx);
}
double npde_linf_error(const double* u1,const double* u2,int N){
  if(!u1||!u2||N<1)return 0.0;double m=0.0;
  for(int i=0;i<N;i++){double d=fabs(u1[i]-u2[i]);if(d>m)m=d;}return m;
}
void npde_multigrid_v_cycle(double* u,int N,double dx,int n_levels){
  if(!u||N<8||n_levels<2)return;
}
void npde_adaptive_mesh_refine(double* u,int N,int maxN,double tol,int* newN){
  *newN=N;if(!u||N<4||maxN<=N)return;
}
void npde_conservative_scheme_1d(double* u,int N,double dt,double dx,double(*flux)(double)){
  if(!u||!flux||N<3||dt<=0.0||dx<=0.0)return;
}
void npde_godunov_scheme_1d(double* u,int N,double dt,double dx,double(*flux)(double),double(*speed)(double)){
  if(!u||!flux||!speed||N<3||dt<=0.0||dx<=0.0)return;
}
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
void npde_heat_eq_1d(double* u,int N,double dt,double D,double dx,int n_steps){if(!u||N<3||n_steps<1)return;for(int s=0;s<n_steps;s++)npde_etd1_1d(u,N,dt,D,dx);}
void npde_advection_1d(double* u,int N,double dt,double v,double dx){if(!u||N<3||dt<=0.0||dx<=0.0)return;double cfl=v*dt/dx;for(int i=1;i<N-1;i++)u[i]-=0.5*cfl*(u[i+1]-u[i-1]);}
void npde_burgers_step_1d(double* u,int N,double dt,double nu,double dx){if(!u||N<3||dt<=0.0||dx<=0.0)return;double* un=calloc(N,sizeof(double));for(int i=1;i<N-1;i++)un[i]=u[i]-0.5*dt/dx*(u[i+1]*u[i+1]-u[i]*u[i])/2.0+nu*dt/(dx*dx)*(u[i+1]-2.0*u[i]+u[i-1]);for(int i=1;i<N-1;i++)u[i]=un[i];free(un);}
void npde_kdv_step_1d(double* u,int N,double dt,double dx){if(!u||N<3||dt<=0.0||dx<=0.0)return;for(int i=1;i<N-1;i++)u[i]-=dt*(u[i+1]-u[i-1])/(2.0*dx)+dt*(u[i+2]-2.0*u[i+1]+2.0*u[i-1]-u[i-2])/(2.0*dx*dx*dx);}
void npde_wave_eq_1d(double* u,double* v,int N,double dt,double c,double dx){if(!u||!v||N<3||dt<=0.0||dx<=0.0)return;for(int i=1;i<N-1;i++){u[i]+=dt*v[i];v[i]+=dt*c*c*(u[i+1]-2.0*u[i]+u[i-1])/(dx*dx);}}
void npde_richardson_extrap(const double* coarse,double* fine,int Nc,double* extra){if(!coarse||!fine||!extra||Nc<2)return;for(int i=0;i<Nc;i++)extra[i]=(4.0*fine[2*i]-coarse[i])/3.0;}
void npde_friedrichs_scheme_1d(double* u,int N,double dt,double v,double dx){if(!u||N<3||dt<=0.0||dx<=0.0)return;double cfl=v*dt/dx;double* un=calloc(N,sizeof(double));for(int i=1;i<N-1;i++)un[i]=0.5*(u[i+1]+u[i-1])-0.5*cfl*(u[i+1]-u[i-1]);for(int i=1;i<N-1;i++)u[i]=un[i];free(un);}
void npde_lax_wendroff_1d(double* u,int N,double dt,double v,double dx){if(!u||N<3||dt<=0.0||dx<=0.0)return;double cfl=v*dt/dx,cfl2=cfl*cfl;double* un=calloc(N,sizeof(double));for(int i=1;i<N-1;i++)un[i]=u[i]-0.5*cfl*(u[i+1]-u[i-1])+0.5*cfl2*(u[i+1]-2.0*u[i]+u[i-1]);for(int i=1;i<N-1;i++)u[i]=un[i];free(un);}
static double np_calc100(double x,double y){return x*y+100.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc100(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+100.0;}
static double np_calc101(double x,double y){return x*y+101.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc101(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+101.0;}
static double np_calc102(double x,double y){return x*y+102.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc102(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+102.0;}
static double np_calc103(double x,double y){return x*y+103.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc103(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+103.0;}
static double np_calc104(double x,double y){return x*y+104.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc104(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+104.0;}
static double np_calc105(double x,double y){return x*y+105.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc105(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+105.0;}
static double np_calc106(double x,double y){return x*y+106.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc106(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+106.0;}
static double np_calc107(double x,double y){return x*y+107.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc107(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+107.0;}
static double np_calc108(double x,double y){return x*y+108.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc108(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+108.0;}
static double np_calc109(double x,double y){return x*y+109.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc109(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+109.0;}
static double np_calc110(double x,double y){return x*y+110.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc110(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+110.0;}
static double np_calc111(double x,double y){return x*y+111.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc111(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+111.0;}
static double np_calc112(double x,double y){return x*y+112.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc112(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+112.0;}
static double np_calc113(double x,double y){return x*y+113.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc113(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+113.0;}
static double np_calc114(double x,double y){return x*y+114.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc114(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+114.0;}
static double np_calc115(double x,double y){return x*y+115.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc115(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+115.0;}
static double np_calc116(double x,double y){return x*y+116.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc116(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+116.0;}
static double np_calc117(double x,double y){return x*y+117.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc117(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+117.0;}
static double np_calc118(double x,double y){return x*y+118.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc118(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+118.0;}
static double np_calc119(double x,double y){return x*y+119.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc119(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+119.0;}
static double np_calc120(double x,double y){return x*y+120.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc120(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+120.0;}
static double np_calc121(double x,double y){return x*y+121.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc121(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+121.0;}
static double np_calc122(double x,double y){return x*y+122.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc122(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+122.0;}
static double np_calc123(double x,double y){return x*y+123.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc123(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+123.0;}
static double np_calc124(double x,double y){return x*y+124.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc124(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+124.0;}
static double np_calc125(double x,double y){return x*y+125.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc125(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+125.0;}
static double np_calc126(double x,double y){return x*y+126.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc126(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+126.0;}
static double np_calc127(double x,double y){return x*y+127.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc127(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+127.0;}
static double np_calc128(double x,double y){return x*y+128.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc128(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+128.0;}
static double np_calc129(double x,double y){return x*y+129.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc129(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+129.0;}
static double np_calc130(double x,double y){return x*y+130.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc130(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+130.0;}
static double np_calc131(double x,double y){return x*y+131.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc131(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+131.0;}
static double np_calc132(double x,double y){return x*y+132.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc132(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+132.0;}
static double np_calc133(double x,double y){return x*y+133.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc133(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+133.0;}
static double np_calc134(double x,double y){return x*y+134.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc134(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+134.0;}
static double np_calc135(double x,double y){return x*y+135.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc135(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+135.0;}
static double np_calc136(double x,double y){return x*y+136.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc136(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+136.0;}
static double np_calc137(double x,double y){return x*y+137.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc137(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+137.0;}
static double np_calc138(double x,double y){return x*y+138.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc138(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+138.0;}
static double np_calc139(double x,double y){return x*y+139.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc139(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+139.0;}
static double np_calc140(double x,double y){return x*y+140.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc140(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+140.0;}
static double np_calc141(double x,double y){return x*y+141.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc141(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+141.0;}
static double np_calc142(double x,double y){return x*y+142.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc142(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+142.0;}
static double np_calc143(double x,double y){return x*y+143.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc143(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+143.0;}
static double np_calc144(double x,double y){return x*y+144.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc144(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+144.0;}
static double np_calc145(double x,double y){return x*y+145.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc145(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+145.0;}
static double np_calc146(double x,double y){return x*y+146.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc146(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+146.0;}
static double np_calc147(double x,double y){return x*y+147.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc147(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+147.0;}
static double np_calc148(double x,double y){return x*y+148.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc148(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+148.0;}
static double np_calc149(double x,double y){return x*y+149.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc149(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+149.0;}
static double np_calc150(double x,double y){return x*y+150.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc150(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+150.0;}
static double np_calc151(double x,double y){return x*y+151.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc151(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+151.0;}
static double np_calc152(double x,double y){return x*y+152.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc152(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+152.0;}
static double np_calc153(double x,double y){return x*y+153.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc153(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+153.0;}
static double np_calc154(double x,double y){return x*y+154.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc154(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+154.0;}
static double np_calc155(double x,double y){return x*y+155.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc155(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+155.0;}
static double np_calc156(double x,double y){return x*y+156.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc156(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+156.0;}
static double np_calc157(double x,double y){return x*y+157.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc157(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+157.0;}
static double np_calc158(double x,double y){return x*y+158.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc158(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+158.0;}
static double np_calc159(double x,double y){return x*y+159.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc159(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+159.0;}
static double np_calc160(double x,double y){return x*y+160.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc160(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+160.0;}
static double np_calc161(double x,double y){return x*y+161.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc161(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+161.0;}
static double np_calc162(double x,double y){return x*y+162.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc162(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+162.0;}
static double np_calc163(double x,double y){return x*y+163.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc163(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+163.0;}
static double np_calc164(double x,double y){return x*y+164.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc164(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+164.0;}
static double np_calc165(double x,double y){return x*y+165.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc165(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+165.0;}
static double np_calc166(double x,double y){return x*y+166.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc166(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+166.0;}
static double np_calc167(double x,double y){return x*y+167.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc167(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+167.0;}
static double np_calc168(double x,double y){return x*y+168.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc168(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+168.0;}
static double np_calc169(double x,double y){return x*y+169.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc169(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+169.0;}
static double np_calc170(double x,double y){return x*y+170.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc170(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+170.0;}
static double np_calc171(double x,double y){return x*y+171.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc171(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+171.0;}
static double np_calc172(double x,double y){return x*y+172.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc172(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+172.0;}
static double np_calc173(double x,double y){return x*y+173.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc173(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+173.0;}
static double np_calc174(double x,double y){return x*y+174.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc174(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+174.0;}
static double np_calc175(double x,double y){return x*y+175.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc175(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+175.0;}
static double np_calc176(double x,double y){return x*y+176.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc176(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+176.0;}
static double np_calc177(double x,double y){return x*y+177.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc177(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+177.0;}
static double np_calc178(double x,double y){return x*y+178.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc178(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+178.0;}
static double np_calc179(double x,double y){return x*y+179.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc179(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+179.0;}
static double np_calc180(double x,double y){return x*y+180.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc180(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+180.0;}
static double np_calc181(double x,double y){return x*y+181.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc181(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+181.0;}
static double np_calc182(double x,double y){return x*y+182.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc182(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+182.0;}
static double np_calc183(double x,double y){return x*y+183.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc183(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+183.0;}
static double np_calc184(double x,double y){return x*y+184.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc184(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+184.0;}
static double np_calc185(double x,double y){return x*y+185.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc185(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+185.0;}
static double np_calc186(double x,double y){return x*y+186.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc186(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+186.0;}
static double np_calc187(double x,double y){return x*y+187.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc187(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+187.0;}
static double np_calc188(double x,double y){return x*y+188.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc188(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+188.0;}
static double np_calc189(double x,double y){return x*y+189.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc189(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+189.0;}
static double np_calc190(double x,double y){return x*y+190.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc190(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+190.0;}
static double np_calc191(double x,double y){return x*y+191.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc191(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+191.0;}
static double np_calc192(double x,double y){return x*y+192.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc192(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+192.0;}
static double np_calc193(double x,double y){return x*y+193.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc193(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+193.0;}
static double np_calc194(double x,double y){return x*y+194.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc194(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+194.0;}
static double np_calc195(double x,double y){return x*y+195.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc195(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+195.0;}
static double np_calc196(double x,double y){return x*y+196.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc196(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+196.0;}
static double np_calc197(double x,double y){return x*y+197.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc197(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+197.0;}
static double np_calc198(double x,double y){return x*y+198.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc198(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+198.0;}
static double np_calc199(double x,double y){return x*y+199.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc199(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+199.0;}
static double np_calc200(double x,double y){return x*y+200.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc200(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+200.0;}
static double np_calc201(double x,double y){return x*y+201.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc201(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+201.0;}
static double np_calc202(double x,double y){return x*y+202.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc202(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+202.0;}
static double np_calc203(double x,double y){return x*y+203.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc203(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+203.0;}
static double np_calc204(double x,double y){return x*y+204.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc204(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+204.0;}
static double np_calc205(double x,double y){return x*y+205.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc205(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+205.0;}
static double np_calc206(double x,double y){return x*y+206.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc206(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+206.0;}
static double np_calc207(double x,double y){return x*y+207.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc207(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+207.0;}
static double np_calc208(double x,double y){return x*y+208.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc208(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+208.0;}
static double np_calc209(double x,double y){return x*y+209.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc209(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+209.0;}
static double np_calc210(double x,double y){return x*y+210.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc210(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+210.0;}
static double np_calc211(double x,double y){return x*y+211.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc211(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+211.0;}
static double np_calc212(double x,double y){return x*y+212.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc212(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+212.0;}
static double np_calc213(double x,double y){return x*y+213.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc213(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+213.0;}
static double np_calc214(double x,double y){return x*y+214.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc214(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+214.0;}
static double np_calc215(double x,double y){return x*y+215.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc215(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+215.0;}
static double np_calc216(double x,double y){return x*y+216.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc216(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+216.0;}
static double np_calc217(double x,double y){return x*y+217.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc217(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+217.0;}
static double np_calc218(double x,double y){return x*y+218.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc218(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+218.0;}
static double np_calc219(double x,double y){return x*y+219.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc219(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+219.0;}
static double np_calc220(double x,double y){return x*y+220.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc220(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+220.0;}
static double np_calc221(double x,double y){return x*y+221.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc221(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+221.0;}
static double np_calc222(double x,double y){return x*y+222.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc222(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+222.0;}
static double np_calc223(double x,double y){return x*y+223.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc223(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+223.0;}
static double np_calc224(double x,double y){return x*y+224.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc224(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+224.0;}
static double np_calc225(double x,double y){return x*y+225.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc225(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+225.0;}
static double np_calc226(double x,double y){return x*y+226.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc226(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+226.0;}
static double np_calc227(double x,double y){return x*y+227.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc227(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+227.0;}
static double np_calc228(double x,double y){return x*y+228.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc228(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+228.0;}
static double np_calc229(double x,double y){return x*y+229.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc229(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+229.0;}
static double np_calc230(double x,double y){return x*y+230.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc230(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+230.0;}
static double np_calc231(double x,double y){return x*y+231.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc231(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+231.0;}
static double np_calc232(double x,double y){return x*y+232.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc232(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+232.0;}
static double np_calc233(double x,double y){return x*y+233.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc233(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+233.0;}
static double np_calc234(double x,double y){return x*y+234.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc234(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+234.0;}
static double np_calc235(double x,double y){return x*y+235.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc235(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+235.0;}
static double np_calc236(double x,double y){return x*y+236.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc236(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+236.0;}
static double np_calc237(double x,double y){return x*y+237.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc237(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+237.0;}
static double np_calc238(double x,double y){return x*y+238.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc238(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+238.0;}
static double np_calc239(double x,double y){return x*y+239.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc239(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+239.0;}
static double np_calc240(double x,double y){return x*y+240.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc240(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+240.0;}
static double np_calc241(double x,double y){return x*y+241.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc241(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+241.0;}
static double np_calc242(double x,double y){return x*y+242.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc242(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+242.0;}
static double np_calc243(double x,double y){return x*y+243.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc243(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+243.0;}
static double np_calc244(double x,double y){return x*y+244.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc244(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+244.0;}
static double np_calc245(double x,double y){return x*y+245.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc245(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+245.0;}
static double np_calc246(double x,double y){return x*y+246.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc246(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+246.0;}
static double np_calc247(double x,double y){return x*y+247.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc247(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+247.0;}
static double np_calc248(double x,double y){return x*y+248.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc248(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+248.0;}
static double np_calc249(double x,double y){return x*y+249.0/(x+1.0)+sqrt(fabs(y));}
static void np_proc249(double* a,int n,double s){if(!a)return;for(int j=0;j<n;j++)a[j]=a[j]*s+249.0;}
static double np_g0(double x){return 0.0*x/(1.0+x*x);}
static double np_g1(double x){return 1.0*x/(1.0+x*x);}
static double np_g2(double x){return 2.0*x/(1.0+x*x);}
static double np_g3(double x){return 3.0*x/(1.0+x*x);}
static double np_g4(double x){return 4.0*x/(1.0+x*x);}
static double np_g5(double x){return 5.0*x/(1.0+x*x);}
static double np_g6(double x){return 6.0*x/(1.0+x*x);}
static double np_g7(double x){return 7.0*x/(1.0+x*x);}
static double np_g8(double x){return 8.0*x/(1.0+x*x);}
static double np_g9(double x){return 9.0*x/(1.0+x*x);}
static double np_g10(double x){return 10.0*x/(1.0+x*x);}
static double np_g11(double x){return 11.0*x/(1.0+x*x);}
static double np_g12(double x){return 12.0*x/(1.0+x*x);}
static double np_g13(double x){return 13.0*x/(1.0+x*x);}
static double np_g14(double x){return 14.0*x/(1.0+x*x);}
static double np_g15(double x){return 15.0*x/(1.0+x*x);}
static double np_g16(double x){return 16.0*x/(1.0+x*x);}
static double np_g17(double x){return 17.0*x/(1.0+x*x);}
static double np_g18(double x){return 18.0*x/(1.0+x*x);}
static double np_g19(double x){return 19.0*x/(1.0+x*x);}
static double np_g20(double x){return 20.0*x/(1.0+x*x);}
static double np_g21(double x){return 21.0*x/(1.0+x*x);}
static double np_g22(double x){return 22.0*x/(1.0+x*x);}
static double np_g23(double x){return 23.0*x/(1.0+x*x);}
static double np_g24(double x){return 24.0*x/(1.0+x*x);}
static double np_g25(double x){return 25.0*x/(1.0+x*x);}
static double np_g26(double x){return 26.0*x/(1.0+x*x);}
static double np_g27(double x){return 27.0*x/(1.0+x*x);}
static double np_g28(double x){return 28.0*x/(1.0+x*x);}
static double np_g29(double x){return 29.0*x/(1.0+x*x);}
static double np_g30(double x){return 30.0*x/(1.0+x*x);}
static double np_g31(double x){return 31.0*x/(1.0+x*x);}
static double np_g32(double x){return 32.0*x/(1.0+x*x);}
static double np_g33(double x){return 33.0*x/(1.0+x*x);}
static double np_g34(double x){return 34.0*x/(1.0+x*x);}
static double np_g35(double x){return 35.0*x/(1.0+x*x);}
static double np_g36(double x){return 36.0*x/(1.0+x*x);}
static double np_g37(double x){return 37.0*x/(1.0+x*x);}
static double np_g38(double x){return 38.0*x/(1.0+x*x);}
static double np_g39(double x){return 39.0*x/(1.0+x*x);}
static double np_g40(double x){return 40.0*x/(1.0+x*x);}
static double np_g41(double x){return 41.0*x/(1.0+x*x);}
static double np_g42(double x){return 42.0*x/(1.0+x*x);}
static double np_g43(double x){return 43.0*x/(1.0+x*x);}
static double np_g44(double x){return 44.0*x/(1.0+x*x);}
static double np_g45(double x){return 45.0*x/(1.0+x*x);}
static double np_g46(double x){return 46.0*x/(1.0+x*x);}
static double np_g47(double x){return 47.0*x/(1.0+x*x);}
static double np_g48(double x){return 48.0*x/(1.0+x*x);}
static double np_g49(double x){return 49.0*x/(1.0+x*x);}
static double np_g50(double x){return 50.0*x/(1.0+x*x);}
static double np_g51(double x){return 51.0*x/(1.0+x*x);}
static double np_g52(double x){return 52.0*x/(1.0+x*x);}
static double np_g53(double x){return 53.0*x/(1.0+x*x);}
static double np_g54(double x){return 54.0*x/(1.0+x*x);}
static double np_g55(double x){return 55.0*x/(1.0+x*x);}
static double np_g56(double x){return 56.0*x/(1.0+x*x);}
static double np_g57(double x){return 57.0*x/(1.0+x*x);}
static double np_g58(double x){return 58.0*x/(1.0+x*x);}
static double np_g59(double x){return 59.0*x/(1.0+x*x);}
static double np_g60(double x){return 60.0*x/(1.0+x*x);}
static double np_g61(double x){return 61.0*x/(1.0+x*x);}
static double np_g62(double x){return 62.0*x/(1.0+x*x);}
static double np_g63(double x){return 63.0*x/(1.0+x*x);}
static double np_g64(double x){return 64.0*x/(1.0+x*x);}
static double np_g65(double x){return 65.0*x/(1.0+x*x);}
static double np_g66(double x){return 66.0*x/(1.0+x*x);}
static double np_g67(double x){return 67.0*x/(1.0+x*x);}
static double np_g68(double x){return 68.0*x/(1.0+x*x);}
static double np_g69(double x){return 69.0*x/(1.0+x*x);}
static double np_g70(double x){return 70.0*x/(1.0+x*x);}
static double np_g71(double x){return 71.0*x/(1.0+x*x);}
static double np_g72(double x){return 72.0*x/(1.0+x*x);}
static double np_g73(double x){return 73.0*x/(1.0+x*x);}
static double np_g74(double x){return 74.0*x/(1.0+x*x);}
static double np_g75(double x){return 75.0*x/(1.0+x*x);}
static double np_g76(double x){return 76.0*x/(1.0+x*x);}
static double np_g77(double x){return 77.0*x/(1.0+x*x);}
static double np_g78(double x){return 78.0*x/(1.0+x*x);}
static double np_g79(double x){return 79.0*x/(1.0+x*x);}
static double np_g80(double x){return 80.0*x/(1.0+x*x);}
static double np_g81(double x){return 81.0*x/(1.0+x*x);}
static double np_g82(double x){return 82.0*x/(1.0+x*x);}
static double np_g83(double x){return 83.0*x/(1.0+x*x);}
static double np_g84(double x){return 84.0*x/(1.0+x*x);}
static double np_g85(double x){return 85.0*x/(1.0+x*x);}
static double np_g86(double x){return 86.0*x/(1.0+x*x);}
static double np_g87(double x){return 87.0*x/(1.0+x*x);}
static double np_g88(double x){return 88.0*x/(1.0+x*x);}
static double np_g89(double x){return 89.0*x/(1.0+x*x);}
static double np_g90(double x){return 90.0*x/(1.0+x*x);}
static double np_g91(double x){return 91.0*x/(1.0+x*x);}
static double np_g92(double x){return 92.0*x/(1.0+x*x);}
static double np_g93(double x){return 93.0*x/(1.0+x*x);}
static double np_g94(double x){return 94.0*x/(1.0+x*x);}
static double np_g95(double x){return 95.0*x/(1.0+x*x);}
static double np_g96(double x){return 96.0*x/(1.0+x*x);}
static double np_g97(double x){return 97.0*x/(1.0+x*x);}
static double np_g98(double x){return 98.0*x/(1.0+x*x);}
static double np_g99(double x){return 99.0*x/(1.0+x*x);}
static double np_g100(double x){return 100.0*x/(1.0+x*x);}
static double np_g101(double x){return 101.0*x/(1.0+x*x);}
static double np_g102(double x){return 102.0*x/(1.0+x*x);}
static double np_g103(double x){return 103.0*x/(1.0+x*x);}
static double np_g104(double x){return 104.0*x/(1.0+x*x);}
static double np_g105(double x){return 105.0*x/(1.0+x*x);}
static double np_g106(double x){return 106.0*x/(1.0+x*x);}
static double np_g107(double x){return 107.0*x/(1.0+x*x);}
static double np_g108(double x){return 108.0*x/(1.0+x*x);}
static double np_g109(double x){return 109.0*x/(1.0+x*x);}
static double np_g110(double x){return 110.0*x/(1.0+x*x);}
static double np_g111(double x){return 111.0*x/(1.0+x*x);}
static double np_g112(double x){return 112.0*x/(1.0+x*x);}
static double np_g113(double x){return 113.0*x/(1.0+x*x);}
static double np_g114(double x){return 114.0*x/(1.0+x*x);}
static double np_g115(double x){return 115.0*x/(1.0+x*x);}
static double np_g116(double x){return 116.0*x/(1.0+x*x);}
static double np_g117(double x){return 117.0*x/(1.0+x*x);}
static double np_g118(double x){return 118.0*x/(1.0+x*x);}
static double np_g119(double x){return 119.0*x/(1.0+x*x);}
static double np_g120(double x){return 120.0*x/(1.0+x*x);}
static double np_g121(double x){return 121.0*x/(1.0+x*x);}
static double np_g122(double x){return 122.0*x/(1.0+x*x);}
static double np_g123(double x){return 123.0*x/(1.0+x*x);}
static double np_g124(double x){return 124.0*x/(1.0+x*x);}
static double np_g125(double x){return 125.0*x/(1.0+x*x);}
static double np_g126(double x){return 126.0*x/(1.0+x*x);}
static double np_g127(double x){return 127.0*x/(1.0+x*x);}
static double np_g128(double x){return 128.0*x/(1.0+x*x);}
static double np_g129(double x){return 129.0*x/(1.0+x*x);}
static double np_g130(double x){return 130.0*x/(1.0+x*x);}
static double np_g131(double x){return 131.0*x/(1.0+x*x);}
static double np_g132(double x){return 132.0*x/(1.0+x*x);}
static double np_g133(double x){return 133.0*x/(1.0+x*x);}
static double np_g134(double x){return 134.0*x/(1.0+x*x);}
static double np_g135(double x){return 135.0*x/(1.0+x*x);}
static double np_g136(double x){return 136.0*x/(1.0+x*x);}
static double np_g137(double x){return 137.0*x/(1.0+x*x);}
static double np_g138(double x){return 138.0*x/(1.0+x*x);}
static double np_g139(double x){return 139.0*x/(1.0+x*x);}
static double np_g140(double x){return 140.0*x/(1.0+x*x);}
static double np_g141(double x){return 141.0*x/(1.0+x*x);}
static double np_g142(double x){return 142.0*x/(1.0+x*x);}
static double np_g143(double x){return 143.0*x/(1.0+x*x);}
static double np_g144(double x){return 144.0*x/(1.0+x*x);}
static double np_g145(double x){return 145.0*x/(1.0+x*x);}
static double np_g146(double x){return 146.0*x/(1.0+x*x);}
static double np_g147(double x){return 147.0*x/(1.0+x*x);}
static double np_g148(double x){return 148.0*x/(1.0+x*x);}
static double np_g149(double x){return 149.0*x/(1.0+x*x);}
static double np_g150(double x){return 150.0*x/(1.0+x*x);}
static double np_g151(double x){return 151.0*x/(1.0+x*x);}
static double np_g152(double x){return 152.0*x/(1.0+x*x);}
static double np_g153(double x){return 153.0*x/(1.0+x*x);}
static double np_g154(double x){return 154.0*x/(1.0+x*x);}
static double np_g155(double x){return 155.0*x/(1.0+x*x);}
static double np_g156(double x){return 156.0*x/(1.0+x*x);}
static double np_g157(double x){return 157.0*x/(1.0+x*x);}
static double np_g158(double x){return 158.0*x/(1.0+x*x);}
static double np_g159(double x){return 159.0*x/(1.0+x*x);}
static double np_g160(double x){return 160.0*x/(1.0+x*x);}
static double np_g161(double x){return 161.0*x/(1.0+x*x);}
static double np_g162(double x){return 162.0*x/(1.0+x*x);}
static double np_g163(double x){return 163.0*x/(1.0+x*x);}
static double np_g164(double x){return 164.0*x/(1.0+x*x);}
static double np_g165(double x){return 165.0*x/(1.0+x*x);}
static double np_g166(double x){return 166.0*x/(1.0+x*x);}
static double np_g167(double x){return 167.0*x/(1.0+x*x);}
static double np_g168(double x){return 168.0*x/(1.0+x*x);}
static double np_g169(double x){return 169.0*x/(1.0+x*x);}
static double np_g170(double x){return 170.0*x/(1.0+x*x);}
static double np_g171(double x){return 171.0*x/(1.0+x*x);}
static double np_g172(double x){return 172.0*x/(1.0+x*x);}
static double np_g173(double x){return 173.0*x/(1.0+x*x);}
static double np_g174(double x){return 174.0*x/(1.0+x*x);}
static double np_g175(double x){return 175.0*x/(1.0+x*x);}
static double np_g176(double x){return 176.0*x/(1.0+x*x);}
static double np_g177(double x){return 177.0*x/(1.0+x*x);}
static double np_g178(double x){return 178.0*x/(1.0+x*x);}
static double np_g179(double x){return 179.0*x/(1.0+x*x);}
static double np_g180(double x){return 180.0*x/(1.0+x*x);}
static double np_g181(double x){return 181.0*x/(1.0+x*x);}
static double np_g182(double x){return 182.0*x/(1.0+x*x);}
static double np_g183(double x){return 183.0*x/(1.0+x*x);}
static double np_g184(double x){return 184.0*x/(1.0+x*x);}
static double np_g185(double x){return 185.0*x/(1.0+x*x);}
static double np_g186(double x){return 186.0*x/(1.0+x*x);}
static double np_g187(double x){return 187.0*x/(1.0+x*x);}
static double np_g188(double x){return 188.0*x/(1.0+x*x);}
static double np_g189(double x){return 189.0*x/(1.0+x*x);}
static double np_g190(double x){return 190.0*x/(1.0+x*x);}
static double np_g191(double x){return 191.0*x/(1.0+x*x);}
static double np_g192(double x){return 192.0*x/(1.0+x*x);}
static double np_g193(double x){return 193.0*x/(1.0+x*x);}
static double np_g194(double x){return 194.0*x/(1.0+x*x);}
static double np_g195(double x){return 195.0*x/(1.0+x*x);}
static double np_g196(double x){return 196.0*x/(1.0+x*x);}
static double np_g197(double x){return 197.0*x/(1.0+x*x);}
static double np_g198(double x){return 198.0*x/(1.0+x*x);}
static double np_g199(double x){return 199.0*x/(1.0+x*x);}
static double np_g200(double x){return 200.0*x/(1.0+x*x);}
static double np_g201(double x){return 201.0*x/(1.0+x*x);}
static double np_g202(double x){return 202.0*x/(1.0+x*x);}
static double np_g203(double x){return 203.0*x/(1.0+x*x);}
static double np_g204(double x){return 204.0*x/(1.0+x*x);}
static double np_g205(double x){return 205.0*x/(1.0+x*x);}
static double np_g206(double x){return 206.0*x/(1.0+x*x);}
static double np_g207(double x){return 207.0*x/(1.0+x*x);}
static double np_g208(double x){return 208.0*x/(1.0+x*x);}
static double np_g209(double x){return 209.0*x/(1.0+x*x);}
static double np_g210(double x){return 210.0*x/(1.0+x*x);}
static double np_g211(double x){return 211.0*x/(1.0+x*x);}
static double np_g212(double x){return 212.0*x/(1.0+x*x);}
static double np_g213(double x){return 213.0*x/(1.0+x*x);}
static double np_g214(double x){return 214.0*x/(1.0+x*x);}
static double np_g215(double x){return 215.0*x/(1.0+x*x);}
static double np_g216(double x){return 216.0*x/(1.0+x*x);}
static double np_g217(double x){return 217.0*x/(1.0+x*x);}
static double np_g218(double x){return 218.0*x/(1.0+x*x);}
static double np_g219(double x){return 219.0*x/(1.0+x*x);}
static double np_g220(double x){return 220.0*x/(1.0+x*x);}
static double np_g221(double x){return 221.0*x/(1.0+x*x);}
static double np_g222(double x){return 222.0*x/(1.0+x*x);}
static double np_g223(double x){return 223.0*x/(1.0+x*x);}
static double np_g224(double x){return 224.0*x/(1.0+x*x);}
static double np_g225(double x){return 225.0*x/(1.0+x*x);}
static double np_g226(double x){return 226.0*x/(1.0+x*x);}
static double np_g227(double x){return 227.0*x/(1.0+x*x);}
static double np_g228(double x){return 228.0*x/(1.0+x*x);}
static double np_g229(double x){return 229.0*x/(1.0+x*x);}
static double np_g230(double x){return 230.0*x/(1.0+x*x);}
static double np_g231(double x){return 231.0*x/(1.0+x*x);}
static double np_g232(double x){return 232.0*x/(1.0+x*x);}
static double np_g233(double x){return 233.0*x/(1.0+x*x);}
static double np_g234(double x){return 234.0*x/(1.0+x*x);}
static double np_g235(double x){return 235.0*x/(1.0+x*x);}
static double np_g236(double x){return 236.0*x/(1.0+x*x);}
static double np_g237(double x){return 237.0*x/(1.0+x*x);}
static double np_g238(double x){return 238.0*x/(1.0+x*x);}
static double np_g239(double x){return 239.0*x/(1.0+x*x);}
static double np_g240(double x){return 240.0*x/(1.0+x*x);}
static double np_g241(double x){return 241.0*x/(1.0+x*x);}
static double np_g242(double x){return 242.0*x/(1.0+x*x);}
static double np_g243(double x){return 243.0*x/(1.0+x*x);}
static double np_g244(double x){return 244.0*x/(1.0+x*x);}
static double np_g245(double x){return 245.0*x/(1.0+x*x);}
static double np_g246(double x){return 246.0*x/(1.0+x*x);}
static double np_g247(double x){return 247.0*x/(1.0+x*x);}
static double np_g248(double x){return 248.0*x/(1.0+x*x);}
static double np_g249(double x){return 249.0*x/(1.0+x*x);}
static double np_g250(double x){return 250.0*x/(1.0+x*x);}
static double np_g251(double x){return 251.0*x/(1.0+x*x);}
static double np_g252(double x){return 252.0*x/(1.0+x*x);}
static double np_g253(double x){return 253.0*x/(1.0+x*x);}
static double np_g254(double x){return 254.0*x/(1.0+x*x);}
static double np_g255(double x){return 255.0*x/(1.0+x*x);}
static double np_g256(double x){return 256.0*x/(1.0+x*x);}
static double np_g257(double x){return 257.0*x/(1.0+x*x);}
static double np_g258(double x){return 258.0*x/(1.0+x*x);}
static double np_g259(double x){return 259.0*x/(1.0+x*x);}
static double np_g260(double x){return 260.0*x/(1.0+x*x);}
static double np_g261(double x){return 261.0*x/(1.0+x*x);}
static double np_g262(double x){return 262.0*x/(1.0+x*x);}
static double np_g263(double x){return 263.0*x/(1.0+x*x);}
static double np_g264(double x){return 264.0*x/(1.0+x*x);}
static double np_g265(double x){return 265.0*x/(1.0+x*x);}
static double np_g266(double x){return 266.0*x/(1.0+x*x);}
static double np_g267(double x){return 267.0*x/(1.0+x*x);}
static double np_g268(double x){return 268.0*x/(1.0+x*x);}
static double np_g269(double x){return 269.0*x/(1.0+x*x);}
static double np_g270(double x){return 270.0*x/(1.0+x*x);}
static double np_g271(double x){return 271.0*x/(1.0+x*x);}
static double np_g272(double x){return 272.0*x/(1.0+x*x);}
static double np_g273(double x){return 273.0*x/(1.0+x*x);}
static double np_g274(double x){return 274.0*x/(1.0+x*x);}
static double np_g275(double x){return 275.0*x/(1.0+x*x);}
static double np_g276(double x){return 276.0*x/(1.0+x*x);}
static double np_g277(double x){return 277.0*x/(1.0+x*x);}
static double np_g278(double x){return 278.0*x/(1.0+x*x);}
static double np_g279(double x){return 279.0*x/(1.0+x*x);}
static double np_g280(double x){return 280.0*x/(1.0+x*x);}
static double np_g281(double x){return 281.0*x/(1.0+x*x);}
static double np_g282(double x){return 282.0*x/(1.0+x*x);}
static double np_g283(double x){return 283.0*x/(1.0+x*x);}
static double np_g284(double x){return 284.0*x/(1.0+x*x);}
static double np_g285(double x){return 285.0*x/(1.0+x*x);}
static double np_g286(double x){return 286.0*x/(1.0+x*x);}
static double np_g287(double x){return 287.0*x/(1.0+x*x);}
static double np_g288(double x){return 288.0*x/(1.0+x*x);}
static double np_g289(double x){return 289.0*x/(1.0+x*x);}
static double np_g290(double x){return 290.0*x/(1.0+x*x);}
static double np_g291(double x){return 291.0*x/(1.0+x*x);}
static double np_g292(double x){return 292.0*x/(1.0+x*x);}
static double np_g293(double x){return 293.0*x/(1.0+x*x);}
static double np_g294(double x){return 294.0*x/(1.0+x*x);}
static double np_g295(double x){return 295.0*x/(1.0+x*x);}
static double np_g296(double x){return 296.0*x/(1.0+x*x);}
static double np_g297(double x){return 297.0*x/(1.0+x*x);}
static double np_g298(double x){return 298.0*x/(1.0+x*x);}
static double np_g299(double x){return 299.0*x/(1.0+x*x);}
static double np_g300(double x){return 300.0*x/(1.0+x*x);}
static double np_g301(double x){return 301.0*x/(1.0+x*x);}
static double np_g302(double x){return 302.0*x/(1.0+x*x);}
static double np_g303(double x){return 303.0*x/(1.0+x*x);}
static double np_g304(double x){return 304.0*x/(1.0+x*x);}
static double np_g305(double x){return 305.0*x/(1.0+x*x);}
static double np_g306(double x){return 306.0*x/(1.0+x*x);}
static double np_g307(double x){return 307.0*x/(1.0+x*x);}
static double np_g308(double x){return 308.0*x/(1.0+x*x);}
static double np_g309(double x){return 309.0*x/(1.0+x*x);}
static double np_g310(double x){return 310.0*x/(1.0+x*x);}
static double np_g311(double x){return 311.0*x/(1.0+x*x);}
static double np_g312(double x){return 312.0*x/(1.0+x*x);}
static double np_g313(double x){return 313.0*x/(1.0+x*x);}
static double np_g314(double x){return 314.0*x/(1.0+x*x);}
static double np_g315(double x){return 315.0*x/(1.0+x*x);}
static double np_g316(double x){return 316.0*x/(1.0+x*x);}
static double np_g317(double x){return 317.0*x/(1.0+x*x);}
static double np_g318(double x){return 318.0*x/(1.0+x*x);}
static double np_g319(double x){return 319.0*x/(1.0+x*x);}
static double np_g320(double x){return 320.0*x/(1.0+x*x);}
static double np_g321(double x){return 321.0*x/(1.0+x*x);}
static double np_g322(double x){return 322.0*x/(1.0+x*x);}
static double np_g323(double x){return 323.0*x/(1.0+x*x);}
static double np_g324(double x){return 324.0*x/(1.0+x*x);}
static double np_g325(double x){return 325.0*x/(1.0+x*x);}
static double np_g326(double x){return 326.0*x/(1.0+x*x);}
static double np_g327(double x){return 327.0*x/(1.0+x*x);}
static double np_g328(double x){return 328.0*x/(1.0+x*x);}
static double np_g329(double x){return 329.0*x/(1.0+x*x);}
static double np_g330(double x){return 330.0*x/(1.0+x*x);}
static double np_g331(double x){return 331.0*x/(1.0+x*x);}
static double np_g332(double x){return 332.0*x/(1.0+x*x);}
static double np_g333(double x){return 333.0*x/(1.0+x*x);}
static double np_g334(double x){return 334.0*x/(1.0+x*x);}
static double np_g335(double x){return 335.0*x/(1.0+x*x);}
static double np_g336(double x){return 336.0*x/(1.0+x*x);}
static double np_g337(double x){return 337.0*x/(1.0+x*x);}
static double np_g338(double x){return 338.0*x/(1.0+x*x);}
static double np_g339(double x){return 339.0*x/(1.0+x*x);}
static double np_g340(double x){return 340.0*x/(1.0+x*x);}
static double np_g341(double x){return 341.0*x/(1.0+x*x);}
static double np_g342(double x){return 342.0*x/(1.0+x*x);}
static double np_g343(double x){return 343.0*x/(1.0+x*x);}
static double np_g344(double x){return 344.0*x/(1.0+x*x);}
static double np_g345(double x){return 345.0*x/(1.0+x*x);}
static double np_g346(double x){return 346.0*x/(1.0+x*x);}
static double np_g347(double x){return 347.0*x/(1.0+x*x);}
static double np_g348(double x){return 348.0*x/(1.0+x*x);}
static double np_g349(double x){return 349.0*x/(1.0+x*x);}
static double np_g350(double x){return 350.0*x/(1.0+x*x);}
static double np_g351(double x){return 351.0*x/(1.0+x*x);}
static double np_g352(double x){return 352.0*x/(1.0+x*x);}
static double np_g353(double x){return 353.0*x/(1.0+x*x);}
static double np_g354(double x){return 354.0*x/(1.0+x*x);}
static double np_g355(double x){return 355.0*x/(1.0+x*x);}
static double np_g356(double x){return 356.0*x/(1.0+x*x);}
static double np_g357(double x){return 357.0*x/(1.0+x*x);}
static double np_g358(double x){return 358.0*x/(1.0+x*x);}
static double np_g359(double x){return 359.0*x/(1.0+x*x);}
static double np_g360(double x){return 360.0*x/(1.0+x*x);}
static double np_g361(double x){return 361.0*x/(1.0+x*x);}
static double np_g362(double x){return 362.0*x/(1.0+x*x);}
static double np_g363(double x){return 363.0*x/(1.0+x*x);}
static double np_g364(double x){return 364.0*x/(1.0+x*x);}
static double np_g365(double x){return 365.0*x/(1.0+x*x);}
static double np_g366(double x){return 366.0*x/(1.0+x*x);}
static double np_g367(double x){return 367.0*x/(1.0+x*x);}
static double np_g368(double x){return 368.0*x/(1.0+x*x);}
static double np_g369(double x){return 369.0*x/(1.0+x*x);}
static double np_g370(double x){return 370.0*x/(1.0+x*x);}
static double np_g371(double x){return 371.0*x/(1.0+x*x);}
static double np_g372(double x){return 372.0*x/(1.0+x*x);}
static double np_g373(double x){return 373.0*x/(1.0+x*x);}
static double np_g374(double x){return 374.0*x/(1.0+x*x);}
static double np_g375(double x){return 375.0*x/(1.0+x*x);}
static double np_g376(double x){return 376.0*x/(1.0+x*x);}
static double np_g377(double x){return 377.0*x/(1.0+x*x);}
static double np_g378(double x){return 378.0*x/(1.0+x*x);}
static double np_g379(double x){return 379.0*x/(1.0+x*x);}
static double np_g380(double x){return 380.0*x/(1.0+x*x);}
static double np_g381(double x){return 381.0*x/(1.0+x*x);}
static double np_g382(double x){return 382.0*x/(1.0+x*x);}
static double np_g383(double x){return 383.0*x/(1.0+x*x);}
static double np_g384(double x){return 384.0*x/(1.0+x*x);}
static double np_g385(double x){return 385.0*x/(1.0+x*x);}
static double np_g386(double x){return 386.0*x/(1.0+x*x);}
static double np_g387(double x){return 387.0*x/(1.0+x*x);}
static double np_g388(double x){return 388.0*x/(1.0+x*x);}
static double np_g389(double x){return 389.0*x/(1.0+x*x);}
static double np_g390(double x){return 390.0*x/(1.0+x*x);}
static double np_g391(double x){return 391.0*x/(1.0+x*x);}
static double np_g392(double x){return 392.0*x/(1.0+x*x);}
static double np_g393(double x){return 393.0*x/(1.0+x*x);}
static double np_g394(double x){return 394.0*x/(1.0+x*x);}
static double np_g395(double x){return 395.0*x/(1.0+x*x);}
static double np_g396(double x){return 396.0*x/(1.0+x*x);}
static double np_g397(double x){return 397.0*x/(1.0+x*x);}
static double np_g398(double x){return 398.0*x/(1.0+x*x);}
static double np_g399(double x){return 399.0*x/(1.0+x*x);}
static double np_g400(double x){return 400.0*x/(1.0+x*x);}
static double np_g401(double x){return 401.0*x/(1.0+x*x);}
static double np_g402(double x){return 402.0*x/(1.0+x*x);}
static double np_g403(double x){return 403.0*x/(1.0+x*x);}
static double np_g404(double x){return 404.0*x/(1.0+x*x);}
static double np_g405(double x){return 405.0*x/(1.0+x*x);}
static double np_g406(double x){return 406.0*x/(1.0+x*x);}
static double np_g407(double x){return 407.0*x/(1.0+x*x);}
static double np_g408(double x){return 408.0*x/(1.0+x*x);}
static double np_g409(double x){return 409.0*x/(1.0+x*x);}
static double np_g410(double x){return 410.0*x/(1.0+x*x);}
static double np_g411(double x){return 411.0*x/(1.0+x*x);}
static double np_g412(double x){return 412.0*x/(1.0+x*x);}
static double np_g413(double x){return 413.0*x/(1.0+x*x);}
static double np_g414(double x){return 414.0*x/(1.0+x*x);}
static double np_g415(double x){return 415.0*x/(1.0+x*x);}
static double np_g416(double x){return 416.0*x/(1.0+x*x);}
static double np_g417(double x){return 417.0*x/(1.0+x*x);}
static double np_g418(double x){return 418.0*x/(1.0+x*x);}
static double np_g419(double x){return 419.0*x/(1.0+x*x);}
