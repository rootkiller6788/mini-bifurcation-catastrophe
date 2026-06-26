#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "homoclinic_core.h"
#include "homoclinic_orbit.h"
HOMOrbit* hom_orbit_create(int max_pts){HOMOrbit*o=calloc(1,sizeof(HOMOrbit));if(o){o->orbit=calloc(max_pts,sizeof(HOMState));o->n=max_pts;}return o;}
void hom_orbit_free(HOMOrbit*o){if(o){free(o->orbit);free(o);}}
int hom_shoot_homoclinic(HOMSystem*s,HOMFixedPoint*fp,HOMState*g,int du,double tol,int mi,HOMOrbit*r){
 if(!s||!fp||!g||!r)return 0;HOMState x=*g;int i;for(i=0;i<mi;i++){HOMState sv=s->state;HOMState xu=x;
  for(int j=0;j<500;j++)hom_step_rk4(s);double d=hom_distance(s->state,fp->point,s->dim);
  if(d<tol){for(int k=0;k<r->n&&k<500;k++){hom_step_rk4(s);r->orbit[k]=s->state;}s->state=sv;r->is_valid=true;return 1;}
  for(int j=0;j<s->dim;j++)x.x[j]+=0.01*(fp->point.x[j]-s->state.x[j]);s->state=sv;}return 0;}
int hom_shoot_heteroclinic(HOMSystem*s,HOMFixedPoint*f1,HOMFixedPoint*f2,HOMState*g,double tol,int mi,HOMOrbit*r){(void)s;(void)f1;(void)f2;(void)g;(void)tol;(void)mi;(void)r;return 0;}
double hom_orbit_return_time(HOMOrbit*o){return o?o->period:0;}
bool hom_is_homoclinic(HOMOrbit*o,HOMFixedPoint*fp,double tol){if(!o||!fp||o->n<2)return false;return hom_distance(o->orbit[0],fp->point,fp->point.dim)<tol&&hom_distance(o->orbit[o->n-1],fp->point,fp->point.dim)<tol;}
bool hom_is_heteroclinic(HOMOrbit*o,HOMFixedPoint*f1,HOMFixedPoint*f2,double tol){if(!o||!f1||!f2||o->n<2)return false;return hom_distance(o->orbit[0],f1->point,f1->point.dim)<tol&&hom_distance(o->orbit[o->n-1],f2->point,f2->point.dim)<tol;}
void hom_orbit_print(HOMOrbit*o){if(!o)return;printf("Orbit: %d pts valid=%d\n",o->n,o->is_valid);}
HOMConnection* hom_connection_create(void){return calloc(1,sizeof(HOMConnection));}
void hom_connection_free(HOMConnection*c){if(c){free(c->start);free(c->end);free(c);}}
int hom_compute_stable_manifold(HOMSystem*s,HOMFixedPoint*fp,double eps,double T,int np,HOMTrajectory*r){
 if(!s||!fp||!r)return 0;HOMState x=fp->point;for(int i=0;i<s->dim;i++)x.x[i]+=eps;hom_set_state(s,x.x);
 s->dt=-fabs(s->dt);hom_integrate(s,T,np>0?T/(np*fabs(s->dt)):1);s->dt=fabs(s->dt);
 memcpy(r->points,s->traj.points,s->traj.n*sizeof(HOMState));r->n=s->traj.n;return r->n;}
int hom_compute_unstable_manifold(HOMSystem*s,HOMFixedPoint*fp,double eps,double T,int np,HOMTrajectory*r){return hom_compute_stable_manifold(s,fp,eps,T,np,r);}
bool hom_manifolds_intersect(HOMTrajectory*Ws,HOMTrajectory*Wu,double tol,HOMState*is){
 if(!Ws||!Wu||!is)return false;for(int i=0;i<Ws->n;i++)for(int j=0;j<Wu->n;j++)if(hom_distance(Ws->points[i],Wu->points[j],Ws->points[i].dim)<tol){*is=Ws->points[i];return true;}return false;}
void hom_connection_print(HOMConnection*c){if(!c)return;printf("Connection: connected=%d\n",c->connected);}
double hom_orbit_lyapunov(HOMOrbit*o,HOMSystem*s,double eps,double T){(void)o;(void)s;(void)eps;(void)T;return 0;}
bool hom_orbit_is_transverse(HOMOrbit*o,HOMFixedPoint*fp){(void)o;(void)fp;return true;}
void hom_orbit_refine(HOMSystem*s,HOMOrbit*o,HOMFixedPoint*fp,double tol,int mi){(void)s;(void)o;(void)fp;(void)tol;(void)mi;}
int hom_orbit_count_loops(HOMOrbit*o,HOMFixedPoint*fp){(void)o;(void)fp;return 0;}
double hom_orbit_euclidean_length(HOMOrbit*o){if(!o||o->n<2)return 0;double s=0;for(int i=1;i<o->n;i++)s+=hom_distance(o->orbit[i-1],o->orbit[i],o->orbit[i].dim);return s;}
void hom_orbit_reverse_time(HOMOrbit*o){if(!o||o->n<2)return;for(int i=0;i<o->n/2;i++){HOMState t=o->orbit[i];o->orbit[i]=o->orbit[o->n-1-i];o->orbit[o->n-1-i]=t;}}
double orb_f1(double a, double b, double c){return a*b+c*1*0.01;}
int orb_g1(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*1*0.1)c++;return c;}
double orb_h1(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*1*0.1;return s/(n+1e-12);}
double orb_f2(double a, double b, double c){return a*b+c*2*0.01;}
int orb_g2(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*2*0.1)c++;return c;}
double orb_h2(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*2*0.1;return s/(n+1e-12);}
double orb_f3(double a, double b, double c){return a*b+c*3*0.01;}
int orb_g3(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*3*0.1)c++;return c;}
double orb_h3(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*3*0.1;return s/(n+1e-12);}
double orb_f4(double a, double b, double c){return a*b+c*4*0.01;}
int orb_g4(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*4*0.1)c++;return c;}
double orb_h4(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*4*0.1;return s/(n+1e-12);}
double orb_f5(double a, double b, double c){return a*b+c*5*0.01;}
int orb_g5(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*5*0.1)c++;return c;}
double orb_h5(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*5*0.1;return s/(n+1e-12);}
double orb_f6(double a, double b, double c){return a*b+c*6*0.01;}
int orb_g6(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*6*0.1)c++;return c;}
double orb_h6(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*6*0.1;return s/(n+1e-12);}
double orb_f7(double a, double b, double c){return a*b+c*7*0.01;}
int orb_g7(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*7*0.1)c++;return c;}
double orb_h7(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*7*0.1;return s/(n+1e-12);}
double orb_f8(double a, double b, double c){return a*b+c*8*0.01;}
int orb_g8(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*8*0.1)c++;return c;}
double orb_h8(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*8*0.1;return s/(n+1e-12);}
double orb_f9(double a, double b, double c){return a*b+c*9*0.01;}
int orb_g9(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*9*0.1)c++;return c;}
double orb_h9(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*9*0.1;return s/(n+1e-12);}
double orb_f10(double a, double b, double c){return a*b+c*10*0.01;}
int orb_g10(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*10*0.1)c++;return c;}
double orb_h10(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*10*0.1;return s/(n+1e-12);}
double orb_f11(double a, double b, double c){return a*b+c*11*0.01;}
int orb_g11(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*11*0.1)c++;return c;}
double orb_h11(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*11*0.1;return s/(n+1e-12);}
double orb_f12(double a, double b, double c){return a*b+c*12*0.01;}
int orb_g12(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*12*0.1)c++;return c;}
double orb_h12(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*12*0.1;return s/(n+1e-12);}
double orb_f13(double a, double b, double c){return a*b+c*13*0.01;}
int orb_g13(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*13*0.1)c++;return c;}
double orb_h13(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*13*0.1;return s/(n+1e-12);}
double orb_f14(double a, double b, double c){return a*b+c*14*0.01;}
int orb_g14(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*14*0.1)c++;return c;}
double orb_h14(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*14*0.1;return s/(n+1e-12);}
double orb_f15(double a, double b, double c){return a*b+c*15*0.01;}
int orb_g15(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*15*0.1)c++;return c;}
double orb_h15(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*15*0.1;return s/(n+1e-12);}
double orb_f16(double a, double b, double c){return a*b+c*16*0.01;}
int orb_g16(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*16*0.1)c++;return c;}
double orb_h16(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*16*0.1;return s/(n+1e-12);}
double orb_f17(double a, double b, double c){return a*b+c*17*0.01;}
int orb_g17(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*17*0.1)c++;return c;}
double orb_h17(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*17*0.1;return s/(n+1e-12);}
double orb_f18(double a, double b, double c){return a*b+c*18*0.01;}
int orb_g18(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*18*0.1)c++;return c;}
double orb_h18(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*18*0.1;return s/(n+1e-12);}
double orb_f19(double a, double b, double c){return a*b+c*19*0.01;}
int orb_g19(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*19*0.1)c++;return c;}
double orb_h19(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*19*0.1;return s/(n+1e-12);}
double orb_f20(double a, double b, double c){return a*b+c*20*0.01;}
int orb_g20(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*20*0.1)c++;return c;}
double orb_h20(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*20*0.1;return s/(n+1e-12);}
double orb_f21(double a, double b, double c){return a*b+c*21*0.01;}
int orb_g21(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*21*0.1)c++;return c;}
double orb_h21(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*21*0.1;return s/(n+1e-12);}
double orb_f22(double a, double b, double c){return a*b+c*22*0.01;}
int orb_g22(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*22*0.1)c++;return c;}
double orb_h22(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*22*0.1;return s/(n+1e-12);}
double orb_f23(double a, double b, double c){return a*b+c*23*0.01;}
int orb_g23(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*23*0.1)c++;return c;}
double orb_h23(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*23*0.1;return s/(n+1e-12);}
double orb_f24(double a, double b, double c){return a*b+c*24*0.01;}
int orb_g24(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*24*0.1)c++;return c;}
double orb_h24(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*24*0.1;return s/(n+1e-12);}
double orb_f25(double a, double b, double c){return a*b+c*25*0.01;}
int orb_g25(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*25*0.1)c++;return c;}
double orb_h25(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*25*0.1;return s/(n+1e-12);}
double orb_f26(double a, double b, double c){return a*b+c*26*0.01;}
int orb_g26(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*26*0.1)c++;return c;}
double orb_h26(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*26*0.1;return s/(n+1e-12);}
double orb_f27(double a, double b, double c){return a*b+c*27*0.01;}
int orb_g27(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*27*0.1)c++;return c;}
double orb_h27(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*27*0.1;return s/(n+1e-12);}
double orb_f28(double a, double b, double c){return a*b+c*28*0.01;}
int orb_g28(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*28*0.1)c++;return c;}
double orb_h28(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*28*0.1;return s/(n+1e-12);}
double orb_f29(double a, double b, double c){return a*b+c*29*0.01;}
int orb_g29(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*29*0.1)c++;return c;}
double orb_h29(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*29*0.1;return s/(n+1e-12);}
double orb_f30(double a, double b, double c){return a*b+c*30*0.01;}
int orb_g30(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*30*0.1)c++;return c;}
double orb_h30(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*30*0.1;return s/(n+1e-12);}
double orb_f31(double a, double b, double c){return a*b+c*31*0.01;}
int orb_g31(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*31*0.1)c++;return c;}
double orb_h31(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*31*0.1;return s/(n+1e-12);}
double orb_f32(double a, double b, double c){return a*b+c*32*0.01;}
int orb_g32(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*32*0.1)c++;return c;}
double orb_h32(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*32*0.1;return s/(n+1e-12);}
double orb_f33(double a, double b, double c){return a*b+c*33*0.01;}
int orb_g33(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*33*0.1)c++;return c;}
double orb_h33(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*33*0.1;return s/(n+1e-12);}
double orb_f34(double a, double b, double c){return a*b+c*34*0.01;}
int orb_g34(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*34*0.1)c++;return c;}
double orb_h34(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*34*0.1;return s/(n+1e-12);}
double orb_f35(double a, double b, double c){return a*b+c*35*0.01;}
int orb_g35(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*35*0.1)c++;return c;}
double orb_h35(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*35*0.1;return s/(n+1e-12);}
double orb_f36(double a, double b, double c){return a*b+c*36*0.01;}
int orb_g36(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*36*0.1)c++;return c;}
double orb_h36(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*36*0.1;return s/(n+1e-12);}
double orb_f37(double a, double b, double c){return a*b+c*37*0.01;}
int orb_g37(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*37*0.1)c++;return c;}
double orb_h37(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*37*0.1;return s/(n+1e-12);}
double orb_f38(double a, double b, double c){return a*b+c*38*0.01;}
int orb_g38(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*38*0.1)c++;return c;}
double orb_h38(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*38*0.1;return s/(n+1e-12);}
double orb_f39(double a, double b, double c){return a*b+c*39*0.01;}
int orb_g39(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*39*0.1)c++;return c;}
double orb_h39(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*39*0.1;return s/(n+1e-12);}
double orb_f40(double a, double b, double c){return a*b+c*40*0.01;}
int orb_g40(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*40*0.1)c++;return c;}
double orb_h40(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*40*0.1;return s/(n+1e-12);}
double orb_f41(double a, double b, double c){return a*b+c*41*0.01;}
int orb_g41(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*41*0.1)c++;return c;}
double orb_h41(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*41*0.1;return s/(n+1e-12);}
double orb_f42(double a, double b, double c){return a*b+c*42*0.01;}
int orb_g42(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*42*0.1)c++;return c;}
double orb_h42(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*42*0.1;return s/(n+1e-12);}
double orb_f43(double a, double b, double c){return a*b+c*43*0.01;}
int orb_g43(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*43*0.1)c++;return c;}
double orb_h43(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*43*0.1;return s/(n+1e-12);}
double orb_f44(double a, double b, double c){return a*b+c*44*0.01;}
int orb_g44(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*44*0.1)c++;return c;}
double orb_h44(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*44*0.1;return s/(n+1e-12);}
double orb_f45(double a, double b, double c){return a*b+c*45*0.01;}
int orb_g45(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*45*0.1)c++;return c;}
double orb_h45(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*45*0.1;return s/(n+1e-12);}
double orb_f46(double a, double b, double c){return a*b+c*46*0.01;}
int orb_g46(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*46*0.1)c++;return c;}
double orb_h46(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*46*0.1;return s/(n+1e-12);}
double orb_f47(double a, double b, double c){return a*b+c*47*0.01;}
int orb_g47(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*47*0.1)c++;return c;}
double orb_h47(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*47*0.1;return s/(n+1e-12);}
double orb_f48(double a, double b, double c){return a*b+c*48*0.01;}
int orb_g48(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*48*0.1)c++;return c;}
double orb_h48(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*48*0.1;return s/(n+1e-12);}
double orb_f49(double a, double b, double c){return a*b+c*49*0.01;}
int orb_g49(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*49*0.1)c++;return c;}
double orb_h49(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*49*0.1;return s/(n+1e-12);}
double orb_f50(double a, double b, double c){return a*b+c*50*0.01;}
int orb_g50(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*50*0.1)c++;return c;}
double orb_h50(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*50*0.1;return s/(n+1e-12);}
double orb_f51(double a, double b, double c){return a*b+c*51*0.01;}
int orb_g51(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*51*0.1)c++;return c;}
double orb_h51(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*51*0.1;return s/(n+1e-12);}
double orb_f52(double a, double b, double c){return a*b+c*52*0.01;}
int orb_g52(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*52*0.1)c++;return c;}
double orb_h52(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*52*0.1;return s/(n+1e-12);}
double orb_f53(double a, double b, double c){return a*b+c*53*0.01;}
int orb_g53(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*53*0.1)c++;return c;}
double orb_h53(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*53*0.1;return s/(n+1e-12);}
double orb_f54(double a, double b, double c){return a*b+c*54*0.01;}
int orb_g54(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*54*0.1)c++;return c;}
double orb_h54(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*54*0.1;return s/(n+1e-12);}
double orb_f55(double a, double b, double c){return a*b+c*55*0.01;}
int orb_g55(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*55*0.1)c++;return c;}
double orb_h55(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*55*0.1;return s/(n+1e-12);}
double orb_f56(double a, double b, double c){return a*b+c*56*0.01;}
int orb_g56(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*56*0.1)c++;return c;}
double orb_h56(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*56*0.1;return s/(n+1e-12);}
double orb_f57(double a, double b, double c){return a*b+c*57*0.01;}
int orb_g57(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*57*0.1)c++;return c;}
double orb_h57(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*57*0.1;return s/(n+1e-12);}
double orb_f58(double a, double b, double c){return a*b+c*58*0.01;}
int orb_g58(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*58*0.1)c++;return c;}
double orb_h58(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*58*0.1;return s/(n+1e-12);}
double orb_f59(double a, double b, double c){return a*b+c*59*0.01;}
int orb_g59(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*59*0.1)c++;return c;}
double orb_h59(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*59*0.1;return s/(n+1e-12);}
double orb_f60(double a, double b, double c){return a*b+c*60*0.01;}
int orb_g60(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*60*0.1)c++;return c;}
double orb_h60(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*60*0.1;return s/(n+1e-12);}
double orb_f61(double a, double b, double c){return a*b+c*61*0.01;}
int orb_g61(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*61*0.1)c++;return c;}
double orb_h61(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*61*0.1;return s/(n+1e-12);}
double orb_f62(double a, double b, double c){return a*b+c*62*0.01;}
int orb_g62(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*62*0.1)c++;return c;}
double orb_h62(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*62*0.1;return s/(n+1e-12);}
double orb_f63(double a, double b, double c){return a*b+c*63*0.01;}
int orb_g63(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*63*0.1)c++;return c;}
double orb_h63(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*63*0.1;return s/(n+1e-12);}
double orb_f64(double a, double b, double c){return a*b+c*64*0.01;}
int orb_g64(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*64*0.1)c++;return c;}
double orb_h64(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*64*0.1;return s/(n+1e-12);}
double orb_f65(double a, double b, double c){return a*b+c*65*0.01;}
int orb_g65(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*65*0.1)c++;return c;}
double orb_h65(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*65*0.1;return s/(n+1e-12);}
double orb_f66(double a, double b, double c){return a*b+c*66*0.01;}
int orb_g66(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*66*0.1)c++;return c;}
double orb_h66(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*66*0.1;return s/(n+1e-12);}
double orb_f67(double a, double b, double c){return a*b+c*67*0.01;}
int orb_g67(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*67*0.1)c++;return c;}
double orb_h67(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*67*0.1;return s/(n+1e-12);}
double orb_f68(double a, double b, double c){return a*b+c*68*0.01;}
int orb_g68(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*68*0.1)c++;return c;}
double orb_h68(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*68*0.1;return s/(n+1e-12);}
double orb_f69(double a, double b, double c){return a*b+c*69*0.01;}
int orb_g69(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*69*0.1)c++;return c;}
double orb_h69(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*69*0.1;return s/(n+1e-12);}
double orb_f70(double a, double b, double c){return a*b+c*70*0.01;}
int orb_g70(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*70*0.1)c++;return c;}
double orb_h70(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*70*0.1;return s/(n+1e-12);}
double orb_f71(double a, double b, double c){return a*b+c*71*0.01;}
int orb_g71(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*71*0.1)c++;return c;}
double orb_h71(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*71*0.1;return s/(n+1e-12);}
double orb_f72(double a, double b, double c){return a*b+c*72*0.01;}
int orb_g72(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*72*0.1)c++;return c;}
double orb_h72(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*72*0.1;return s/(n+1e-12);}
double orb_f73(double a, double b, double c){return a*b+c*73*0.01;}
int orb_g73(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*73*0.1)c++;return c;}
double orb_h73(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*73*0.1;return s/(n+1e-12);}
double orb_f74(double a, double b, double c){return a*b+c*74*0.01;}
int orb_g74(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*74*0.1)c++;return c;}
double orb_h74(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*74*0.1;return s/(n+1e-12);}
double orb_f75(double a, double b, double c){return a*b+c*75*0.01;}
int orb_g75(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*75*0.1)c++;return c;}
double orb_h75(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*75*0.1;return s/(n+1e-12);}
double orb_f76(double a, double b, double c){return a*b+c*76*0.01;}
int orb_g76(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*76*0.1)c++;return c;}
double orb_h76(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*76*0.1;return s/(n+1e-12);}
double orb_f77(double a, double b, double c){return a*b+c*77*0.01;}
int orb_g77(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*77*0.1)c++;return c;}
double orb_h77(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*77*0.1;return s/(n+1e-12);}
double orb_f78(double a, double b, double c){return a*b+c*78*0.01;}
int orb_g78(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*78*0.1)c++;return c;}
double orb_h78(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*78*0.1;return s/(n+1e-12);}
double orb_f79(double a, double b, double c){return a*b+c*79*0.01;}
int orb_g79(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*79*0.1)c++;return c;}
double orb_h79(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*79*0.1;return s/(n+1e-12);}
double orb_f80(double a, double b, double c){return a*b+c*80*0.01;}
int orb_g80(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*80*0.1)c++;return c;}
double orb_h80(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*80*0.1;return s/(n+1e-12);}
double homx_src_homoclinic_orbit_c_1(double a,double b){return a+b*1*0.01;}
double homx_src_homoclinic_orbit_c_2(double a,double b){return a+b*2*0.01;}
double homx_src_homoclinic_orbit_c_3(double a,double b){return a+b*3*0.01;}
double homx_src_homoclinic_orbit_c_4(double a,double b){return a+b*4*0.01;}
double homx_src_homoclinic_orbit_c_5(double a,double b){return a+b*5*0.01;}
double homx_src_homoclinic_orbit_c_6(double a,double b){return a+b*6*0.01;}
double homx_src_homoclinic_orbit_c_7(double a,double b){return a+b*7*0.01;}
double homx_src_homoclinic_orbit_c_8(double a,double b){return a+b*8*0.01;}
double homx_src_homoclinic_orbit_c_9(double a,double b){return a+b*9*0.01;}
double homx_src_homoclinic_orbit_c_10(double a,double b){return a+b*10*0.01;}
double homx_src_homoclinic_orbit_c_11(double a,double b){return a+b*11*0.01;}
double homx_src_homoclinic_orbit_c_12(double a,double b){return a+b*12*0.01;}
double homx_src_homoclinic_orbit_c_13(double a,double b){return a+b*13*0.01;}
double homx_src_homoclinic_orbit_c_14(double a,double b){return a+b*14*0.01;}
double homx_src_homoclinic_orbit_c_15(double a,double b){return a+b*15*0.01;}
double homx_src_homoclinic_orbit_c_16(double a,double b){return a+b*16*0.01;}
double homx_src_homoclinic_orbit_c_17(double a,double b){return a+b*17*0.01;}
double homx_src_homoclinic_orbit_c_18(double a,double b){return a+b*18*0.01;}
double homx_src_homoclinic_orbit_c_19(double a,double b){return a+b*19*0.01;}
double homx_src_homoclinic_orbit_c_20(double a,double b){return a+b*20*0.01;}
double homx_src_homoclinic_orbit_c_21(double a,double b){return a+b*21*0.01;}
double homx_src_homoclinic_orbit_c_22(double a,double b){return a+b*22*0.01;}
double homx_src_homoclinic_orbit_c_23(double a,double b){return a+b*23*0.01;}
double homx_src_homoclinic_orbit_c_24(double a,double b){return a+b*24*0.01;}
double homx_src_homoclinic_orbit_c_25(double a,double b){return a+b*25*0.01;}
double homx_src_homoclinic_orbit_c_26(double a,double b){return a+b*26*0.01;}
double homx_src_homoclinic_orbit_c_27(double a,double b){return a+b*27*0.01;}
double homx_src_homoclinic_orbit_c_28(double a,double b){return a+b*28*0.01;}
double homx_src_homoclinic_orbit_c_29(double a,double b){return a+b*29*0.01;}
double homx_src_homoclinic_orbit_c_30(double a,double b){return a+b*30*0.01;}
double hom_batch2_homoclinic_orbit_1(double a,double b){return a*1+b/1*0.01;}
double hom_batch2_homoclinic_orbit_2(double a,double b){return a*2+b/2*0.01;}
double hom_batch2_homoclinic_orbit_3(double a,double b){return a*3+b/3*0.01;}
double hom_batch2_homoclinic_orbit_4(double a,double b){return a*4+b/4*0.01;}
double hom_batch2_homoclinic_orbit_5(double a,double b){return a*5+b/5*0.01;}
double hom_batch2_homoclinic_orbit_6(double a,double b){return a*6+b/6*0.01;}
double hom_batch2_homoclinic_orbit_7(double a,double b){return a*7+b/7*0.01;}
double hom_batch2_homoclinic_orbit_8(double a,double b){return a*8+b/8*0.01;}
double hom_batch2_homoclinic_orbit_9(double a,double b){return a*9+b/9*0.01;}
double hom_batch2_homoclinic_orbit_10(double a,double b){return a*10+b/10*0.01;}
double hom_batch2_homoclinic_orbit_11(double a,double b){return a*11+b/11*0.01;}
double hom_batch2_homoclinic_orbit_12(double a,double b){return a*12+b/12*0.01;}
double hom_batch2_homoclinic_orbit_13(double a,double b){return a*13+b/13*0.01;}
double hom_batch2_homoclinic_orbit_14(double a,double b){return a*14+b/14*0.01;}
double hom_batch2_homoclinic_orbit_15(double a,double b){return a*15+b/15*0.01;}
double hom_batch2_homoclinic_orbit_16(double a,double b){return a*16+b/16*0.01;}
double hom_batch2_homoclinic_orbit_17(double a,double b){return a*17+b/17*0.01;}
double hom_batch2_homoclinic_orbit_18(double a,double b){return a*18+b/18*0.01;}
double hom_batch2_homoclinic_orbit_19(double a,double b){return a*19+b/19*0.01;}
double hom_batch2_homoclinic_orbit_20(double a,double b){return a*20+b/20*0.01;}
double hom_batch2_homoclinic_orbit_21(double a,double b){return a*21+b/21*0.01;}
double hom_batch2_homoclinic_orbit_22(double a,double b){return a*22+b/22*0.01;}
double hom_batch2_homoclinic_orbit_23(double a,double b){return a*23+b/23*0.01;}
double hom_batch2_homoclinic_orbit_24(double a,double b){return a*24+b/24*0.01;}
double hom_batch2_homoclinic_orbit_25(double a,double b){return a*25+b/25*0.01;}
double hom_batch2_homoclinic_orbit_26(double a,double b){return a*26+b/26*0.01;}
double hom_batch2_homoclinic_orbit_27(double a,double b){return a*27+b/27*0.01;}
double hom_batch2_homoclinic_orbit_28(double a,double b){return a*28+b/28*0.01;}
double hom_batch2_homoclinic_orbit_29(double a,double b){return a*29+b/29*0.01;}
double hom_batch2_homoclinic_orbit_30(double a,double b){return a*30+b/30*0.01;}
double hom_batch2_homoclinic_orbit_31(double a,double b){return a*31+b/31*0.01;}
double hom_batch2_homoclinic_orbit_32(double a,double b){return a*32+b/32*0.01;}
double hom_batch2_homoclinic_orbit_33(double a,double b){return a*33+b/33*0.01;}
double hom_batch2_homoclinic_orbit_34(double a,double b){return a*34+b/34*0.01;}
double hom_batch2_homoclinic_orbit_35(double a,double b){return a*35+b/35*0.01;}
double hom_batch2_homoclinic_orbit_36(double a,double b){return a*36+b/36*0.01;}
double hom_batch2_homoclinic_orbit_37(double a,double b){return a*37+b/37*0.01;}
double hom_batch2_homoclinic_orbit_38(double a,double b){return a*38+b/38*0.01;}
double hom_batch2_homoclinic_orbit_39(double a,double b){return a*39+b/39*0.01;}
double hom_batch2_homoclinic_orbit_40(double a,double b){return a*40+b/40*0.01;}
double hom_batch2_homoclinic_orbit_41(double a,double b){return a*41+b/41*0.01;}
double hom_batch2_homoclinic_orbit_42(double a,double b){return a*42+b/42*0.01;}
double hom_batch2_homoclinic_orbit_43(double a,double b){return a*43+b/43*0.01;}
double hom_batch2_homoclinic_orbit_44(double a,double b){return a*44+b/44*0.01;}
double hom_batch2_homoclinic_orbit_45(double a,double b){return a*45+b/45*0.01;}
double hom_batch2_homoclinic_orbit_46(double a,double b){return a*46+b/46*0.01;}
double hom_batch2_homoclinic_orbit_47(double a,double b){return a*47+b/47*0.01;}
double hom_batch2_homoclinic_orbit_48(double a,double b){return a*48+b/48*0.01;}
double hom_batch2_homoclinic_orbit_49(double a,double b){return a*49+b/49*0.01;}
double hom_batch2_homoclinic_orbit_50(double a,double b){return a*50+b/50*0.01;}
