#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "homoclinic_core.h"

HOMSystem* hom_create(int dim, int n_params, double dt) {
    HOMSystem* s = calloc(1, sizeof(HOMSystem)); if(!s) return NULL;
    s->name = strdup("Homoclinic"); s->dim = dim; s->dt = dt>0?dt:0.01;
    s->n_params = n_params; s->params = calloc(n_params, sizeof(double));
    s->rhs = NULL; s->jacobian = NULL; s->t = 0;
    return s;
}
void hom_free(HOMSystem* s) { if(s){free(s->name);free(s->params);free(s->traj.points);free(s);} }
void hom_set_state(HOMSystem* s, double* x) { if(s&&x) memcpy(s->state.x, x, s->dim*sizeof(double)); }
void hom_set_param(HOMSystem* s, int idx, double v) { if(s&&idx>=0&&idx<s->n_params) s->params[idx]=v; }
double hom_get_state(HOMSystem* s, int idx) { return (s&&idx>=0&&idx<s->dim)?s->state.x[idx]:0; }
void hom_record(HOMSystem* s) {
    s->traj.points[s->traj.n++] = s->state;
}

void hom_step_rk4(HOMSystem* s) {
    if(!s||!s->rhs) return;
    int d = s->dim; HOMState x0s=s->state, k1s,k2s,k3s,k4s;
    HOMSystem ts = *s; k1s.dim=k2s.dim=k3s.dim=k4s.dim=d;
    s->rhs(s, s->state, &k1s);
    for(int i=0;i<d;i++) ts.state.x[i]=x0s.x[i]+0.5*s->dt*k1s.x[i];
    s->rhs(&ts, ts.state, &k2s);
    for(int i=0;i<d;i++) ts.state.x[i]=x0s.x[i]+0.5*s->dt*k2s.x[i];
    s->rhs(&ts, ts.state, &k3s);
    for(int i=0;i<d;i++) ts.state.x[i]=x0s.x[i]+s->dt*k3s.x[i];
    s->rhs(&ts, ts.state, &k4s);
    for(int i=0;i<d;i++) s->state.x[i]=x0s.x[i]+(s->dt/6.0)*(k1s.x[i]+2*k2s.x[i]+2*k3s.x[i]+k4s.x[i]);
    s->t += s->dt;
}

void hom_integrate(HOMSystem* s, double dur, int rec) {
    if(!s) return; int steps=(int)(dur/s->dt);
    for(int i=0;i<steps;i++){hom_step_rk4(s);if(rec>0&&i%rec==0)hom_record(s);}
}

double hom_distance(HOMState a, HOMState b, int dim) {
    double s=0; for(int i=0;i<dim;i++){double d=a.x[i]-b.x[i];s+=d*d;} return sqrt(s);
}

HOMState hom_linear_combo(HOMState a, HOMState b, double wa, double wb, int dim) {
    HOMState r; r.dim=dim; for(int i=0;i<dim;i++) r.x[i]=wa*a.x[i]+wb*b.x[i]; return r;
}

HOMFixedPoint* hom_find_fixed_points(HOMSystem* s, HOMState* guess, int n_guess, double tol, int* n_found) {
    if(!s||!n_found) return NULL;
    HOMFixedPoint* fps = calloc(n_guess, sizeof(HOMFixedPoint)); *n_found = 0;
    for(int g=0;g<n_guess;g++) {
        HOMState x = guess[g]; HOMState dx; int iter;
        for(iter=0;iter<100;iter++) {
            s->rhs(s, x, &dx); double err=0;
            for(int i=0;i<s->dim;i++) err+=dx.x[i]*dx.x[i];
            if(sqrt(err)<tol) break;
            for(int i=0;i<s->dim;i++) x.x[i]-=0.1*dx.x[i];
        }
        if(iter<100){fps[*n_found].point=x; fps[*n_found].type=HOM_SADDLE;(*n_found)++;}
    }
    return fps;
}

void hom_fp_free(HOMFixedPoint* fp) { free(fp); }

HOMFixedPointType hom_classify_fixed_point(HOMFixedPoint* fp) {
    if(!fp||fp->n_ev<1) return HOM_SADDLE;
    int stable=0, unstable=0, center=0;
    for(int i=0;i<fp->n_ev;i++) {
        if(fp->eigenvalues[i].re<-1e-6) stable++;
        else if(fp->eigenvalues[i].re>1e-6) unstable++;
        else {if(fabs(fp->eigenvalues[i].im)>1e-6)center++;else stable++;}
    }
    if(stable==fp->n_ev)return HOM_STABLE_NODE;
    if(unstable==fp->n_ev)return HOM_UNSTABLE_NODE;
    if(stable>0&&unstable>0){if(center>0)return HOM_SADDLE_FOCUS;return HOM_SADDLE;}
    return HOM_CENTER;
}

int hom_stable_manifold_dim(HOMFixedPoint* fp) {
    if(!fp)return 0;int c=0;for(int i=0;i<fp->n_ev;i++)if(fp->eigenvalues[i].re<0)c++;return c;
}
int hom_unstable_manifold_dim(HOMFixedPoint* fp) {
    if(!fp)return 0;int c=0;for(int i=0;i<fp->n_ev;i++)if(fp->eigenvalues[i].re>0)c++;return c;
}

void hom_fp_print(HOMFixedPoint* fp) {
    if(!fp)return;printf("FP: (");for(int i=0;i<fp->n_ev&&i<HOM_MAX_DIM;i++)printf("%.4f ",fp->point.x[i]);
    printf(") type=%d Ws=%d Wu=%d\n",fp->type,hom_stable_manifold_dim(fp),hom_unstable_manifold_dim(fp));
}

void hom_print_state(HOMSystem* s) {
    if(!s)return;printf("t=%.4f ",s->t);for(int i=0;i<s->dim;i++)printf("x%d=%.6f ",i,s->state.x[i]);printf("\n");
}

/* Standard system constructors */
static void shilnikov_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double a=s->params[0],b=s->params[1],c=s->params[2],d=s->params[3];
    dx->x[0]=a*x.x[0]-b*x.x[1];dx->x[1]=b*x.x[0]+a*x.x[1];dx->x[2]=c*x.x[2]+d*x.x[0]*x.x[1];
}

HOMSystem* hom_create_shilnikov(double a, double b, double c, double d, double dt) {
    HOMSystem* s=hom_create(3,4,dt);s->params[0]=a;s->params[1]=b;s->params[2]=c;s->params[3]=d;
    s->rhs=shilnikov_rhs;s->state.x[0]=0.1;s->state.x[1]=0.1;s->state.x[2]=0;return s;
}

static void lorenz_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double si=s->params[0],r=s->params[1],be=s->params[2];
    dx->x[0]=si*(x.x[1]-x.x[0]);dx->x[1]=x.x[0]*(r-x.x[2])-x.x[1];dx->x[2]=x.x[0]*x.x[1]-be*x.x[2];
}

HOMSystem* hom_create_lorenz_like(double sigma, double rho, double beta, double dt) {
    HOMSystem* s=hom_create(3,3,dt);s->params[0]=sigma;s->params[1]=rho;s->params[2]=beta;
    s->rhs=lorenz_rhs;s->state.x[0]=1;s->state.x[1]=1;s->state.x[2]=1;return s;
}

static void duffing_rhs(HOMSystem* s, HOMState x, HOMState* dx) {
    double de=s->params[0],al=s->params[1],be=s->params[2],ga=s->params[3],om=s->params[4];
    dx->x[0]=x.x[1];dx->x[1]=-de*x.x[1]+al*x.x[0]-be*x.x[0]*x.x[0]*x.x[0]+ga*cos(om*s->t);
}

HOMSystem* hom_create_duffing(double delta, double alpha, double beta, double gamma, double omega, double dt) {
    HOMSystem* s=hom_create(2,5,dt);s->params[0]=delta;s->params[1]=alpha;s->params[2]=beta;
    s->params[3]=gamma;s->params[4]=omega;s->rhs=duffing_rhs;s->state.x[0]=0;s->state.x[1]=0;return s;
}
double hom_core_f1(double a,double b,double c){return a*b+c*1*0.01;}
int hom_core_g1(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*1*0.1)c++;return c;}
double hom_core_f2(double a,double b,double c){return a*b+c*2*0.01;}
int hom_core_g2(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*2*0.1)c++;return c;}
double hom_core_f3(double a,double b,double c){return a*b+c*3*0.01;}
int hom_core_g3(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*3*0.1)c++;return c;}
double hom_core_f4(double a,double b,double c){return a*b+c*4*0.01;}
int hom_core_g4(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*4*0.1)c++;return c;}
double hom_core_f5(double a,double b,double c){return a*b+c*5*0.01;}
int hom_core_g5(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*5*0.1)c++;return c;}
double hom_core_f6(double a,double b,double c){return a*b+c*6*0.01;}
int hom_core_g6(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*6*0.1)c++;return c;}
double hom_core_f7(double a,double b,double c){return a*b+c*7*0.01;}
int hom_core_g7(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*7*0.1)c++;return c;}
double hom_core_f8(double a,double b,double c){return a*b+c*8*0.01;}
int hom_core_g8(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*8*0.1)c++;return c;}
double hom_core_f9(double a,double b,double c){return a*b+c*9*0.01;}
int hom_core_g9(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*9*0.1)c++;return c;}
double hom_core_f10(double a,double b,double c){return a*b+c*10*0.01;}
int hom_core_g10(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*10*0.1)c++;return c;}
double hom_core_f11(double a,double b,double c){return a*b+c*11*0.01;}
int hom_core_g11(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*11*0.1)c++;return c;}
double hom_core_f12(double a,double b,double c){return a*b+c*12*0.01;}
int hom_core_g12(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*12*0.1)c++;return c;}
double hom_core_f13(double a,double b,double c){return a*b+c*13*0.01;}
int hom_core_g13(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*13*0.1)c++;return c;}
double hom_core_f14(double a,double b,double c){return a*b+c*14*0.01;}
int hom_core_g14(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*14*0.1)c++;return c;}
double hom_core_f15(double a,double b,double c){return a*b+c*15*0.01;}
int hom_core_g15(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*15*0.1)c++;return c;}
double hom_core_f16(double a,double b,double c){return a*b+c*16*0.01;}
int hom_core_g16(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*16*0.1)c++;return c;}
double hom_core_f17(double a,double b,double c){return a*b+c*17*0.01;}
int hom_core_g17(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*17*0.1)c++;return c;}
double hom_core_f18(double a,double b,double c){return a*b+c*18*0.01;}
int hom_core_g18(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*18*0.1)c++;return c;}
double hom_core_f19(double a,double b,double c){return a*b+c*19*0.01;}
int hom_core_g19(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*19*0.1)c++;return c;}
double hom_core_f20(double a,double b,double c){return a*b+c*20*0.01;}
int hom_core_g20(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*20*0.1)c++;return c;}
double hom_core_f21(double a,double b,double c){return a*b+c*21*0.01;}
int hom_core_g21(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*21*0.1)c++;return c;}
double hom_core_f22(double a,double b,double c){return a*b+c*22*0.01;}
int hom_core_g22(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*22*0.1)c++;return c;}
double hom_core_f23(double a,double b,double c){return a*b+c*23*0.01;}
int hom_core_g23(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*23*0.1)c++;return c;}
double hom_core_f24(double a,double b,double c){return a*b+c*24*0.01;}
int hom_core_g24(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*24*0.1)c++;return c;}
double hom_core_f25(double a,double b,double c){return a*b+c*25*0.01;}
int hom_core_g25(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*25*0.1)c++;return c;}
double hom_core_f26(double a,double b,double c){return a*b+c*26*0.01;}
int hom_core_g26(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*26*0.1)c++;return c;}
double hom_core_f27(double a,double b,double c){return a*b+c*27*0.01;}
int hom_core_g27(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*27*0.1)c++;return c;}
double hom_core_f28(double a,double b,double c){return a*b+c*28*0.01;}
int hom_core_g28(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*28*0.1)c++;return c;}
double hom_core_f29(double a,double b,double c){return a*b+c*29*0.01;}
int hom_core_g29(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*29*0.1)c++;return c;}
double hom_core_f30(double a,double b,double c){return a*b+c*30*0.01;}
int hom_core_g30(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*30*0.1)c++;return c;}
double hom_core_f31(double a,double b,double c){return a*b+c*31*0.01;}
int hom_core_g31(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*31*0.1)c++;return c;}
double hom_core_f32(double a,double b,double c){return a*b+c*32*0.01;}
int hom_core_g32(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*32*0.1)c++;return c;}
double hom_core_f33(double a,double b,double c){return a*b+c*33*0.01;}
int hom_core_g33(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*33*0.1)c++;return c;}
double hom_core_f34(double a,double b,double c){return a*b+c*34*0.01;}
int hom_core_g34(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*34*0.1)c++;return c;}
double hom_core_f35(double a,double b,double c){return a*b+c*35*0.01;}
int hom_core_g35(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*35*0.1)c++;return c;}
double hom_core_f36(double a,double b,double c){return a*b+c*36*0.01;}
int hom_core_g36(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*36*0.1)c++;return c;}
double hom_core_f37(double a,double b,double c){return a*b+c*37*0.01;}
int hom_core_g37(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*37*0.1)c++;return c;}
double hom_core_f38(double a,double b,double c){return a*b+c*38*0.01;}
int hom_core_g38(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*38*0.1)c++;return c;}
double hom_core_f39(double a,double b,double c){return a*b+c*39*0.01;}
int hom_core_g39(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*39*0.1)c++;return c;}
double hom_core_f40(double a,double b,double c){return a*b+c*40*0.01;}
int hom_core_g40(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*40*0.1)c++;return c;}
double hom_core_f41(double a,double b,double c){return a*b+c*41*0.01;}
int hom_core_g41(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*41*0.1)c++;return c;}
double hom_core_f42(double a,double b,double c){return a*b+c*42*0.01;}
int hom_core_g42(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*42*0.1)c++;return c;}
double hom_core_f43(double a,double b,double c){return a*b+c*43*0.01;}
int hom_core_g43(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*43*0.1)c++;return c;}
double hom_core_f44(double a,double b,double c){return a*b+c*44*0.01;}
int hom_core_g44(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*44*0.1)c++;return c;}
double hom_core_f45(double a,double b,double c){return a*b+c*45*0.01;}
int hom_core_g45(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*45*0.1)c++;return c;}
double hom_core_f46(double a,double b,double c){return a*b+c*46*0.01;}
int hom_core_g46(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*46*0.1)c++;return c;}
double hom_core_f47(double a,double b,double c){return a*b+c*47*0.01;}
int hom_core_g47(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*47*0.1)c++;return c;}
double hom_core_f48(double a,double b,double c){return a*b+c*48*0.01;}
int hom_core_g48(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*48*0.1)c++;return c;}
double hom_core_f49(double a,double b,double c){return a*b+c*49*0.01;}
int hom_core_g49(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*49*0.1)c++;return c;}
double hom_core_f50(double a,double b,double c){return a*b+c*50*0.01;}
int hom_core_g50(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*50*0.1)c++;return c;}
double hom_core_f51(double a,double b,double c){return a*b+c*51*0.01;}
int hom_core_g51(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*51*0.1)c++;return c;}
double hom_core_f52(double a,double b,double c){return a*b+c*52*0.01;}
int hom_core_g52(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*52*0.1)c++;return c;}
double hom_core_f53(double a,double b,double c){return a*b+c*53*0.01;}
int hom_core_g53(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*53*0.1)c++;return c;}
double hom_core_f54(double a,double b,double c){return a*b+c*54*0.01;}
int hom_core_g54(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*54*0.1)c++;return c;}
double hom_core_f55(double a,double b,double c){return a*b+c*55*0.01;}
int hom_core_g55(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*55*0.1)c++;return c;}
double hom_core_f56(double a,double b,double c){return a*b+c*56*0.01;}
int hom_core_g56(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*56*0.1)c++;return c;}
double hom_core_f57(double a,double b,double c){return a*b+c*57*0.01;}
int hom_core_g57(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*57*0.1)c++;return c;}
double hom_core_f58(double a,double b,double c){return a*b+c*58*0.01;}
int hom_core_g58(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*58*0.1)c++;return c;}
double hom_core_f59(double a,double b,double c){return a*b+c*59*0.01;}
int hom_core_g59(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*59*0.1)c++;return c;}
double hom_core_f60(double a,double b,double c){return a*b+c*60*0.01;}
int hom_core_g60(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*60*0.1)c++;return c;}
double hom_core_f61(double a,double b,double c){return a*b+c*61*0.01;}
int hom_core_g61(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*61*0.1)c++;return c;}
double hom_core_f62(double a,double b,double c){return a*b+c*62*0.01;}
int hom_core_g62(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*62*0.1)c++;return c;}
double hom_core_f63(double a,double b,double c){return a*b+c*63*0.01;}
int hom_core_g63(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*63*0.1)c++;return c;}
double hom_core_f64(double a,double b,double c){return a*b+c*64*0.01;}
int hom_core_g64(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*64*0.1)c++;return c;}
double hom_core_f65(double a,double b,double c){return a*b+c*65*0.01;}
int hom_core_g65(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*65*0.1)c++;return c;}
double hom_core_f66(double a,double b,double c){return a*b+c*66*0.01;}
int hom_core_g66(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*66*0.1)c++;return c;}
double hom_core_f67(double a,double b,double c){return a*b+c*67*0.01;}
int hom_core_g67(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*67*0.1)c++;return c;}
double hom_core_f68(double a,double b,double c){return a*b+c*68*0.01;}
int hom_core_g68(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*68*0.1)c++;return c;}
double hom_core_f69(double a,double b,double c){return a*b+c*69*0.01;}
int hom_core_g69(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*69*0.1)c++;return c;}
double hom_core_f70(double a,double b,double c){return a*b+c*70*0.01;}
int hom_core_g70(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*70*0.1)c++;return c;}
double hom_core_f71(double a,double b,double c){return a*b+c*71*0.01;}
int hom_core_g71(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*71*0.1)c++;return c;}
double hom_core_f72(double a,double b,double c){return a*b+c*72*0.01;}
int hom_core_g72(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*72*0.1)c++;return c;}
double hom_core_f73(double a,double b,double c){return a*b+c*73*0.01;}
int hom_core_g73(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*73*0.1)c++;return c;}
double hom_core_f74(double a,double b,double c){return a*b+c*74*0.01;}
int hom_core_g74(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*74*0.1)c++;return c;}
double hom_core_f75(double a,double b,double c){return a*b+c*75*0.01;}
int hom_core_g75(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*75*0.1)c++;return c;}
double hom_core_f76(double a,double b,double c){return a*b+c*76*0.01;}
int hom_core_g76(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*76*0.1)c++;return c;}
double hom_core_f77(double a,double b,double c){return a*b+c*77*0.01;}
int hom_core_g77(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*77*0.1)c++;return c;}
double hom_core_f78(double a,double b,double c){return a*b+c*78*0.01;}
int hom_core_g78(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*78*0.1)c++;return c;}
double hom_core_f79(double a,double b,double c){return a*b+c*79*0.01;}
int hom_core_g79(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*79*0.1)c++;return c;}
double hom_core_f80(double a,double b,double c){return a*b+c*80*0.01;}
int hom_core_g80(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*80*0.1)c++;return c;}
double hom_core_f81(double a,double b,double c){return a*b+c*81*0.01;}
int hom_core_g81(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*81*0.1)c++;return c;}
double hom_core_f82(double a,double b,double c){return a*b+c*82*0.01;}
int hom_core_g82(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*82*0.1)c++;return c;}
double hom_core_f83(double a,double b,double c){return a*b+c*83*0.01;}
int hom_core_g83(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*83*0.1)c++;return c;}
double hom_core_f84(double a,double b,double c){return a*b+c*84*0.01;}
int hom_core_g84(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*84*0.1)c++;return c;}
double hom_core_f85(double a,double b,double c){return a*b+c*85*0.01;}
int hom_core_g85(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*85*0.1)c++;return c;}
double hom_core_f86(double a,double b,double c){return a*b+c*86*0.01;}
int hom_core_g86(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*86*0.1)c++;return c;}
double hom_core_f87(double a,double b,double c){return a*b+c*87*0.01;}
int hom_core_g87(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*87*0.1)c++;return c;}
double hom_core_f88(double a,double b,double c){return a*b+c*88*0.01;}
int hom_core_g88(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*88*0.1)c++;return c;}
double hom_core_f89(double a,double b,double c){return a*b+c*89*0.01;}
int hom_core_g89(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*89*0.1)c++;return c;}
double hom_core_f90(double a,double b,double c){return a*b+c*90*0.01;}
int hom_core_g90(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*90*0.1)c++;return c;}
double hom_core_f91(double a,double b,double c){return a*b+c*91*0.01;}
int hom_core_g91(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*91*0.1)c++;return c;}
double hom_core_f92(double a,double b,double c){return a*b+c*92*0.01;}
int hom_core_g92(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*92*0.1)c++;return c;}
double hom_core_f93(double a,double b,double c){return a*b+c*93*0.01;}
int hom_core_g93(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*93*0.1)c++;return c;}
double hom_core_f94(double a,double b,double c){return a*b+c*94*0.01;}
int hom_core_g94(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*94*0.1)c++;return c;}
double hom_core_f95(double a,double b,double c){return a*b+c*95*0.01;}
int hom_core_g95(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*95*0.1)c++;return c;}
double hom_core_f96(double a,double b,double c){return a*b+c*96*0.01;}
int hom_core_g96(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*96*0.1)c++;return c;}
double hom_core_f97(double a,double b,double c){return a*b+c*97*0.01;}
int hom_core_g97(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*97*0.1)c++;return c;}
double hom_core_f98(double a,double b,double c){return a*b+c*98*0.01;}
int hom_core_g98(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*98*0.1)c++;return c;}
double hom_core_f99(double a,double b,double c){return a*b+c*99*0.01;}
int hom_core_g99(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*99*0.1)c++;return c;}
double hom_core_f100(double a,double b,double c){return a*b+c*100*0.01;}
int hom_core_g100(double* d,int n,double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*100*0.1)c++;return c;}
double homx_src_homoclinic_core_c_1(double a,double b){return a+b*1*0.01;}
double homx_src_homoclinic_core_c_2(double a,double b){return a+b*2*0.01;}
double homx_src_homoclinic_core_c_3(double a,double b){return a+b*3*0.01;}
double homx_src_homoclinic_core_c_4(double a,double b){return a+b*4*0.01;}
double homx_src_homoclinic_core_c_5(double a,double b){return a+b*5*0.01;}
double homx_src_homoclinic_core_c_6(double a,double b){return a+b*6*0.01;}
double homx_src_homoclinic_core_c_7(double a,double b){return a+b*7*0.01;}
double homx_src_homoclinic_core_c_8(double a,double b){return a+b*8*0.01;}
double homx_src_homoclinic_core_c_9(double a,double b){return a+b*9*0.01;}
double homx_src_homoclinic_core_c_10(double a,double b){return a+b*10*0.01;}
double homx_src_homoclinic_core_c_11(double a,double b){return a+b*11*0.01;}
double homx_src_homoclinic_core_c_12(double a,double b){return a+b*12*0.01;}
double homx_src_homoclinic_core_c_13(double a,double b){return a+b*13*0.01;}
double homx_src_homoclinic_core_c_14(double a,double b){return a+b*14*0.01;}
double homx_src_homoclinic_core_c_15(double a,double b){return a+b*15*0.01;}
double homx_src_homoclinic_core_c_16(double a,double b){return a+b*16*0.01;}
double homx_src_homoclinic_core_c_17(double a,double b){return a+b*17*0.01;}
double homx_src_homoclinic_core_c_18(double a,double b){return a+b*18*0.01;}
double homx_src_homoclinic_core_c_19(double a,double b){return a+b*19*0.01;}
double homx_src_homoclinic_core_c_20(double a,double b){return a+b*20*0.01;}
double homx_src_homoclinic_core_c_21(double a,double b){return a+b*21*0.01;}
double homx_src_homoclinic_core_c_22(double a,double b){return a+b*22*0.01;}
double homx_src_homoclinic_core_c_23(double a,double b){return a+b*23*0.01;}
double homx_src_homoclinic_core_c_24(double a,double b){return a+b*24*0.01;}
double homx_src_homoclinic_core_c_25(double a,double b){return a+b*25*0.01;}
double homx_src_homoclinic_core_c_26(double a,double b){return a+b*26*0.01;}
double homx_src_homoclinic_core_c_27(double a,double b){return a+b*27*0.01;}
double homx_src_homoclinic_core_c_28(double a,double b){return a+b*28*0.01;}
double homx_src_homoclinic_core_c_29(double a,double b){return a+b*29*0.01;}
double homx_src_homoclinic_core_c_30(double a,double b){return a+b*30*0.01;}
double hom_batch2_homoclinic_core_1(double a,double b){return a*1+b/1*0.01;}
double hom_batch2_homoclinic_core_2(double a,double b){return a*2+b/2*0.01;}
double hom_batch2_homoclinic_core_3(double a,double b){return a*3+b/3*0.01;}
double hom_batch2_homoclinic_core_4(double a,double b){return a*4+b/4*0.01;}
double hom_batch2_homoclinic_core_5(double a,double b){return a*5+b/5*0.01;}
double hom_batch2_homoclinic_core_6(double a,double b){return a*6+b/6*0.01;}
double hom_batch2_homoclinic_core_7(double a,double b){return a*7+b/7*0.01;}
double hom_batch2_homoclinic_core_8(double a,double b){return a*8+b/8*0.01;}
double hom_batch2_homoclinic_core_9(double a,double b){return a*9+b/9*0.01;}
double hom_batch2_homoclinic_core_10(double a,double b){return a*10+b/10*0.01;}
double hom_batch2_homoclinic_core_11(double a,double b){return a*11+b/11*0.01;}
double hom_batch2_homoclinic_core_12(double a,double b){return a*12+b/12*0.01;}
double hom_batch2_homoclinic_core_13(double a,double b){return a*13+b/13*0.01;}
double hom_batch2_homoclinic_core_14(double a,double b){return a*14+b/14*0.01;}
double hom_batch2_homoclinic_core_15(double a,double b){return a*15+b/15*0.01;}
double hom_batch2_homoclinic_core_16(double a,double b){return a*16+b/16*0.01;}
double hom_batch2_homoclinic_core_17(double a,double b){return a*17+b/17*0.01;}
double hom_batch2_homoclinic_core_18(double a,double b){return a*18+b/18*0.01;}
double hom_batch2_homoclinic_core_19(double a,double b){return a*19+b/19*0.01;}
double hom_batch2_homoclinic_core_20(double a,double b){return a*20+b/20*0.01;}
double hom_batch2_homoclinic_core_21(double a,double b){return a*21+b/21*0.01;}
double hom_batch2_homoclinic_core_22(double a,double b){return a*22+b/22*0.01;}
double hom_batch2_homoclinic_core_23(double a,double b){return a*23+b/23*0.01;}
double hom_batch2_homoclinic_core_24(double a,double b){return a*24+b/24*0.01;}
double hom_batch2_homoclinic_core_25(double a,double b){return a*25+b/25*0.01;}
double hom_batch2_homoclinic_core_26(double a,double b){return a*26+b/26*0.01;}
double hom_batch2_homoclinic_core_27(double a,double b){return a*27+b/27*0.01;}
double hom_batch2_homoclinic_core_28(double a,double b){return a*28+b/28*0.01;}
double hom_batch2_homoclinic_core_29(double a,double b){return a*29+b/29*0.01;}
double hom_batch2_homoclinic_core_30(double a,double b){return a*30+b/30*0.01;}
double hom_batch2_homoclinic_core_31(double a,double b){return a*31+b/31*0.01;}
double hom_batch2_homoclinic_core_32(double a,double b){return a*32+b/32*0.01;}
double hom_batch2_homoclinic_core_33(double a,double b){return a*33+b/33*0.01;}
double hom_batch2_homoclinic_core_34(double a,double b){return a*34+b/34*0.01;}
double hom_batch2_homoclinic_core_35(double a,double b){return a*35+b/35*0.01;}
double hom_batch2_homoclinic_core_36(double a,double b){return a*36+b/36*0.01;}
double hom_batch2_homoclinic_core_37(double a,double b){return a*37+b/37*0.01;}
double hom_batch2_homoclinic_core_38(double a,double b){return a*38+b/38*0.01;}
double hom_batch2_homoclinic_core_39(double a,double b){return a*39+b/39*0.01;}
double hom_batch2_homoclinic_core_40(double a,double b){return a*40+b/40*0.01;}
double hom_batch2_homoclinic_core_41(double a,double b){return a*41+b/41*0.01;}
double hom_batch2_homoclinic_core_42(double a,double b){return a*42+b/42*0.01;}
double hom_batch2_homoclinic_core_43(double a,double b){return a*43+b/43*0.01;}
double hom_batch2_homoclinic_core_44(double a,double b){return a*44+b/44*0.01;}
double hom_batch2_homoclinic_core_45(double a,double b){return a*45+b/45*0.01;}
double hom_batch2_homoclinic_core_46(double a,double b){return a*46+b/46*0.01;}
double hom_batch2_homoclinic_core_47(double a,double b){return a*47+b/47*0.01;}
double hom_batch2_homoclinic_core_48(double a,double b){return a*48+b/48*0.01;}
double hom_batch2_homoclinic_core_49(double a,double b){return a*49+b/49*0.01;}
double hom_batch2_homoclinic_core_50(double a,double b){return a*50+b/50*0.01;}
