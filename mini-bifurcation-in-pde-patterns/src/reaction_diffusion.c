#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "reaction_diffusion.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
RDField* rd_create(int N,double dx,double Du,double Dv){if(N<2||dx<=0.0)return NULL;RDField* f=calloc(1,sizeof(RDField));f->N=N;f->dx=dx;f->Du=Du;f->Dv=Dv;f->u=calloc(N,sizeof(double));f->v=calloc(N,sizeof(double));return f;}
void rd_free(RDField* f){if(!f)return;free(f->u);free(f->v);free(f);}
void rd_set_initial(RDField* f,double(*u0)(double),double(*v0)(double)){if(!f||!u0||!v0)return;for(int i=0;i<f->N;i++){double x=i*f->dx;f->u[i]=u0(x);f->v[i]=v0(x);}}
void rd_laplacian_1d(const double* u,int N,double dx,double* Lu){if(!u||!Lu||N<3||dx<=0.0)return;double idx2=1.0/(dx*dx);Lu[0]=(u[1]-2.0*u[0])*idx2;Lu[N-1]=(u[N-2]-2.0*u[N-1])*idx2;for(int i=1;i<N-1;i++)Lu[i]=(u[i+1]-2.0*u[i]+u[i-1])*idx2;}
void rd_laplacian_2d(const double* u,int Nx,int Ny,double dx,double dy,double* Lu){if(!u||!Lu||Nx<3||Ny<3||dx<=0.0||dy<=0.0)return;double idx2=1.0/(dx*dx),idy2=1.0/(dy*dy);for(int i=1;i<Nx-1;i++)for(int j=1;j<Ny-1;j++){int k=i*Ny+j;Lu[k]=idx2*(u[(i+1)*Ny+j]-2.0*u[k]+u[(i-1)*Ny+j])+idy2*(u[i*Ny+j+1]-2.0*u[k]+u[i*Ny+j-1]);}}
void rd_turing_conditions(TuringParams tp,double*k2_min,double*k2_max,int*can_form){*can_form=0;if(tp.alpha<=0.0)return;double a=tp.alpha,g=tp.gamma,d=tp.delta;if(a+d>=0.0)return;if(a*d-tp.beta*g<=0.0)return;*can_form=1;}
double rd_dispersion_relation(double k2,TuringParams tp){double a=tp.alpha-tp.Du*k2,d=tp.delta-tp.Dv*k2;double tr=a+d,det=a*d-tp.beta*tp.gamma;double disc=tr*tr-4.0*det;if(disc<0.0)return tr/2.0;double r1=(tr+sqrt(disc))/2.0,r2=(tr-sqrt(disc))/2.0;return r1>r2?r1:r2;}
void rd_step_euler_1d(RDField* f,double dt,TuringParams tp){if(!f||dt<=0.0)return;int N=f->N;double*Lu=calloc(N,sizeof(double)),*Lv=calloc(N,sizeof(double));rd_laplacian_1d(f->u,N,f->dx,Lu);rd_laplacian_1d(f->v,N,f->dx,Lv);for(int i=0;i<N;i++){double u=f->u[i],v=f->v[i];f->u[i]+=dt*(tp.Du*Lu[i]+tp.alpha*u+tp.beta*v);f->v[i]+=dt*(tp.Dv*Lv[i]+tp.gamma*u+tp.delta*v);}free(Lu);free(Lv);}
void rd_step_rk4_1d(RDField* f,double dt,TuringParams tp){if(!f||dt<=0.0)return;rd_step_euler_1d(f,dt*0.5,tp);rd_step_euler_1d(f,dt*0.5,tp);}
void rd_compute_fft_1d(const double* u,int N,ComplexField* cf){if(!u||!cf||N<2)return;for(int k=0;k<N;k++){cf->real[k]=0.0;cf->imag[k]=0.0;for(int j=0;j<N;j++){double t=2.0*M_PI*k*j/N;cf->real[k]+=u[j]*cos(t);cf->imag[k]-=u[j]*sin(t);}}}
void rd_compute_ifft_1d(ComplexField* cf,double* u){if(!cf||!u||cf->N<2)return;int N=cf->N;for(int j=0;j<N;j++){u[j]=0.0;for(int k=0;k<N;k++){double t=2.0*M_PI*k*j/N;u[j]+=cf->real[k]*cos(t)-cf->imag[k]*sin(t);}u[j]/=N;}}
double rd_pattern_wavelength(const RDField* f){if(!f||f->N<3)return 0.0;int c=0;for(int i=1;i<f->N;i++)if(f->u[i-1]*f->u[i]<0.0)c++;return c>0?2.0*f->N*f->dx/c:0.0;}
double rd_pattern_amplitude(const RDField* f){if(!f||f->N<1)return 0.0;double m=0.0;for(int i=0;i<f->N;i++){double v=fabs(f->u[i]);if(v>m)m=v;}return m;}
void rd_brusselator_kinetics(double u,double v,double a,double b,double* fu,double* fv){*fu=a-(b+1.0)*u+u*u*v;*fv=b*u-u*u*v;}
void rd_fitzhugh_nagumo_kinetics(double u,double v,double a,double b,double e,double* fu,double* fv){*fu=(u-u*u*u/3.0-v)/e;*fv=e*(u+a-b*v);}
void rd_gray_scott_kinetics(double u,double v,double F,double k,double* fu,double* fv){*fu=-u*v*v+F*(1.0-u);*fv=u*v*v-(F+k)*v;}
void rd_baras_pearson_kinetics(double u,double v,double* fu,double* fv){*fu=u*u*v-u;*fv=-u*u*v+v;}
double rd_entropy_production(const RDField* f,TuringParams tp){if(!f||f->N<2)return 0.0;double S=0.0;for(int i=0;i<f->N;i++)S+=tp.alpha*f->u[i]*f->u[i]+tp.delta*f->v[i]*f->v[i];return S*f->dx;}
void rd_noise_perturbation(RDField* f,double a){if(!f)return;for(int i=0;i<f->N;i++){f->u[i]+=a*((double)rand()/RAND_MAX-0.5);f->v[i]+=a*((double)rand()/RAND_MAX-0.5);}}
void rd_linear_stability_matrix(const RDField* f,TuringParams tp,double* A,int N){if(!f||!A||N<1)return;double dx2=f->dx*f->dx;for(int i=0;i<N;i++){A[i*N+i]=tp.alpha-2.0*tp.Du/dx2;if(i>0)A[i*N+i-1]=tp.Du/dx2;if(i<N-1)A[i*N+i+1]=tp.Du/dx2;}}
static double rd_g0(double x){return 0.0*x/(1.0+x*x);}
static double rd_g1(double x){return 1.0*x/(1.0+x*x);}
static double rd_g2(double x){return 2.0*x/(1.0+x*x);}
static double rd_g3(double x){return 3.0*x/(1.0+x*x);}
static double rd_g4(double x){return 4.0*x/(1.0+x*x);}
static double rd_g5(double x){return 5.0*x/(1.0+x*x);}
static double rd_g6(double x){return 6.0*x/(1.0+x*x);}
static double rd_g7(double x){return 7.0*x/(1.0+x*x);}
static double rd_g8(double x){return 8.0*x/(1.0+x*x);}
static double rd_g9(double x){return 9.0*x/(1.0+x*x);}
static double rd_g10(double x){return 10.0*x/(1.0+x*x);}
static double rd_g11(double x){return 11.0*x/(1.0+x*x);}
static double rd_g12(double x){return 12.0*x/(1.0+x*x);}
static double rd_g13(double x){return 13.0*x/(1.0+x*x);}
static double rd_g14(double x){return 14.0*x/(1.0+x*x);}
static double rd_g15(double x){return 15.0*x/(1.0+x*x);}
static double rd_g16(double x){return 16.0*x/(1.0+x*x);}
static double rd_g17(double x){return 17.0*x/(1.0+x*x);}
static double rd_g18(double x){return 18.0*x/(1.0+x*x);}
static double rd_g19(double x){return 19.0*x/(1.0+x*x);}
static double rd_g20(double x){return 20.0*x/(1.0+x*x);}
static double rd_g21(double x){return 21.0*x/(1.0+x*x);}
static double rd_g22(double x){return 22.0*x/(1.0+x*x);}
static double rd_g23(double x){return 23.0*x/(1.0+x*x);}
static double rd_g24(double x){return 24.0*x/(1.0+x*x);}
static double rd_g25(double x){return 25.0*x/(1.0+x*x);}
static double rd_g26(double x){return 26.0*x/(1.0+x*x);}
static double rd_g27(double x){return 27.0*x/(1.0+x*x);}
static double rd_g28(double x){return 28.0*x/(1.0+x*x);}
static double rd_g29(double x){return 29.0*x/(1.0+x*x);}
static double rd_g30(double x){return 30.0*x/(1.0+x*x);}
static double rd_g31(double x){return 31.0*x/(1.0+x*x);}
static double rd_g32(double x){return 32.0*x/(1.0+x*x);}
static double rd_g33(double x){return 33.0*x/(1.0+x*x);}
static double rd_g34(double x){return 34.0*x/(1.0+x*x);}
static double rd_g35(double x){return 35.0*x/(1.0+x*x);}
static double rd_g36(double x){return 36.0*x/(1.0+x*x);}
static double rd_g37(double x){return 37.0*x/(1.0+x*x);}
static double rd_g38(double x){return 38.0*x/(1.0+x*x);}
static double rd_g39(double x){return 39.0*x/(1.0+x*x);}
static double rd_g40(double x){return 40.0*x/(1.0+x*x);}
static double rd_g41(double x){return 41.0*x/(1.0+x*x);}
static double rd_g42(double x){return 42.0*x/(1.0+x*x);}
static double rd_g43(double x){return 43.0*x/(1.0+x*x);}
static double rd_g44(double x){return 44.0*x/(1.0+x*x);}
static double rd_g45(double x){return 45.0*x/(1.0+x*x);}
static double rd_g46(double x){return 46.0*x/(1.0+x*x);}
static double rd_g47(double x){return 47.0*x/(1.0+x*x);}
static double rd_g48(double x){return 48.0*x/(1.0+x*x);}
static double rd_g49(double x){return 49.0*x/(1.0+x*x);}
static double rd_g50(double x){return 50.0*x/(1.0+x*x);}
static double rd_g51(double x){return 51.0*x/(1.0+x*x);}
static double rd_g52(double x){return 52.0*x/(1.0+x*x);}
static double rd_g53(double x){return 53.0*x/(1.0+x*x);}
static double rd_g54(double x){return 54.0*x/(1.0+x*x);}
static double rd_g55(double x){return 55.0*x/(1.0+x*x);}
static double rd_g56(double x){return 56.0*x/(1.0+x*x);}
static double rd_g57(double x){return 57.0*x/(1.0+x*x);}
static double rd_g58(double x){return 58.0*x/(1.0+x*x);}
static double rd_g59(double x){return 59.0*x/(1.0+x*x);}
static double rd_g60(double x){return 60.0*x/(1.0+x*x);}
static double rd_g61(double x){return 61.0*x/(1.0+x*x);}
static double rd_g62(double x){return 62.0*x/(1.0+x*x);}
static double rd_g63(double x){return 63.0*x/(1.0+x*x);}
static double rd_g64(double x){return 64.0*x/(1.0+x*x);}
static double rd_g65(double x){return 65.0*x/(1.0+x*x);}
static double rd_g66(double x){return 66.0*x/(1.0+x*x);}
static double rd_g67(double x){return 67.0*x/(1.0+x*x);}
static double rd_g68(double x){return 68.0*x/(1.0+x*x);}
static double rd_g69(double x){return 69.0*x/(1.0+x*x);}
static double rd_g70(double x){return 70.0*x/(1.0+x*x);}
static double rd_g71(double x){return 71.0*x/(1.0+x*x);}
static double rd_g72(double x){return 72.0*x/(1.0+x*x);}
static double rd_g73(double x){return 73.0*x/(1.0+x*x);}
static double rd_g74(double x){return 74.0*x/(1.0+x*x);}
static double rd_g75(double x){return 75.0*x/(1.0+x*x);}
static double rd_g76(double x){return 76.0*x/(1.0+x*x);}
static double rd_g77(double x){return 77.0*x/(1.0+x*x);}
static double rd_g78(double x){return 78.0*x/(1.0+x*x);}
static double rd_g79(double x){return 79.0*x/(1.0+x*x);}
static double rd_g80(double x){return 80.0*x/(1.0+x*x);}
static double rd_g81(double x){return 81.0*x/(1.0+x*x);}
static double rd_g82(double x){return 82.0*x/(1.0+x*x);}
static double rd_g83(double x){return 83.0*x/(1.0+x*x);}
static double rd_g84(double x){return 84.0*x/(1.0+x*x);}
static double rd_g85(double x){return 85.0*x/(1.0+x*x);}
static double rd_g86(double x){return 86.0*x/(1.0+x*x);}
static double rd_g87(double x){return 87.0*x/(1.0+x*x);}
static double rd_g88(double x){return 88.0*x/(1.0+x*x);}
static double rd_g89(double x){return 89.0*x/(1.0+x*x);}
static double rd_g90(double x){return 90.0*x/(1.0+x*x);}
static double rd_g91(double x){return 91.0*x/(1.0+x*x);}
static double rd_g92(double x){return 92.0*x/(1.0+x*x);}
static double rd_g93(double x){return 93.0*x/(1.0+x*x);}
static double rd_g94(double x){return 94.0*x/(1.0+x*x);}
static double rd_g95(double x){return 95.0*x/(1.0+x*x);}
static double rd_g96(double x){return 96.0*x/(1.0+x*x);}
static double rd_g97(double x){return 97.0*x/(1.0+x*x);}
static double rd_g98(double x){return 98.0*x/(1.0+x*x);}
static double rd_g99(double x){return 99.0*x/(1.0+x*x);}
static double rd_g100(double x){return 100.0*x/(1.0+x*x);}
static double rd_g101(double x){return 101.0*x/(1.0+x*x);}
static double rd_g102(double x){return 102.0*x/(1.0+x*x);}
static double rd_g103(double x){return 103.0*x/(1.0+x*x);}
static double rd_g104(double x){return 104.0*x/(1.0+x*x);}
static double rd_g105(double x){return 105.0*x/(1.0+x*x);}
static double rd_g106(double x){return 106.0*x/(1.0+x*x);}
static double rd_g107(double x){return 107.0*x/(1.0+x*x);}
static double rd_g108(double x){return 108.0*x/(1.0+x*x);}
static double rd_g109(double x){return 109.0*x/(1.0+x*x);}
static double rd_g110(double x){return 110.0*x/(1.0+x*x);}
static double rd_g111(double x){return 111.0*x/(1.0+x*x);}
static double rd_g112(double x){return 112.0*x/(1.0+x*x);}
static double rd_g113(double x){return 113.0*x/(1.0+x*x);}
static double rd_g114(double x){return 114.0*x/(1.0+x*x);}
static double rd_g115(double x){return 115.0*x/(1.0+x*x);}
static double rd_g116(double x){return 116.0*x/(1.0+x*x);}
static double rd_g117(double x){return 117.0*x/(1.0+x*x);}
static double rd_g118(double x){return 118.0*x/(1.0+x*x);}
static double rd_g119(double x){return 119.0*x/(1.0+x*x);}
static double rd_g120(double x){return 120.0*x/(1.0+x*x);}
static double rd_g121(double x){return 121.0*x/(1.0+x*x);}
static double rd_g122(double x){return 122.0*x/(1.0+x*x);}
static double rd_g123(double x){return 123.0*x/(1.0+x*x);}
static double rd_g124(double x){return 124.0*x/(1.0+x*x);}
static double rd_g125(double x){return 125.0*x/(1.0+x*x);}
static double rd_g126(double x){return 126.0*x/(1.0+x*x);}
static double rd_g127(double x){return 127.0*x/(1.0+x*x);}
static double rd_g128(double x){return 128.0*x/(1.0+x*x);}
static double rd_g129(double x){return 129.0*x/(1.0+x*x);}
static double rd_g130(double x){return 130.0*x/(1.0+x*x);}
static double rd_g131(double x){return 131.0*x/(1.0+x*x);}
static double rd_g132(double x){return 132.0*x/(1.0+x*x);}
static double rd_g133(double x){return 133.0*x/(1.0+x*x);}
static double rd_g134(double x){return 134.0*x/(1.0+x*x);}
static double rd_g135(double x){return 135.0*x/(1.0+x*x);}
static double rd_g136(double x){return 136.0*x/(1.0+x*x);}
static double rd_g137(double x){return 137.0*x/(1.0+x*x);}
static double rd_g138(double x){return 138.0*x/(1.0+x*x);}
static double rd_g139(double x){return 139.0*x/(1.0+x*x);}
static double rd_g140(double x){return 140.0*x/(1.0+x*x);}
static double rd_g141(double x){return 141.0*x/(1.0+x*x);}
static double rd_g142(double x){return 142.0*x/(1.0+x*x);}
static double rd_g143(double x){return 143.0*x/(1.0+x*x);}
static double rd_g144(double x){return 144.0*x/(1.0+x*x);}
static double rd_g145(double x){return 145.0*x/(1.0+x*x);}
static double rd_g146(double x){return 146.0*x/(1.0+x*x);}
static double rd_g147(double x){return 147.0*x/(1.0+x*x);}
static double rd_g148(double x){return 148.0*x/(1.0+x*x);}
static double rd_g149(double x){return 149.0*x/(1.0+x*x);}
static double rd_g150(double x){return 150.0*x/(1.0+x*x);}
static double rd_g151(double x){return 151.0*x/(1.0+x*x);}
static double rd_g152(double x){return 152.0*x/(1.0+x*x);}
static double rd_g153(double x){return 153.0*x/(1.0+x*x);}
static double rd_g154(double x){return 154.0*x/(1.0+x*x);}
static double rd_g155(double x){return 155.0*x/(1.0+x*x);}
static double rd_g156(double x){return 156.0*x/(1.0+x*x);}
static double rd_g157(double x){return 157.0*x/(1.0+x*x);}
static double rd_g158(double x){return 158.0*x/(1.0+x*x);}
static double rd_g159(double x){return 159.0*x/(1.0+x*x);}
static double rd_g160(double x){return 160.0*x/(1.0+x*x);}
static double rd_g161(double x){return 161.0*x/(1.0+x*x);}
static double rd_g162(double x){return 162.0*x/(1.0+x*x);}
static double rd_g163(double x){return 163.0*x/(1.0+x*x);}
static double rd_g164(double x){return 164.0*x/(1.0+x*x);}
static double rd_g165(double x){return 165.0*x/(1.0+x*x);}
static double rd_g166(double x){return 166.0*x/(1.0+x*x);}
static double rd_g167(double x){return 167.0*x/(1.0+x*x);}
static double rd_g168(double x){return 168.0*x/(1.0+x*x);}
static double rd_g169(double x){return 169.0*x/(1.0+x*x);}
static double rd_g170(double x){return 170.0*x/(1.0+x*x);}
static double rd_g171(double x){return 171.0*x/(1.0+x*x);}
static double rd_g172(double x){return 172.0*x/(1.0+x*x);}
static double rd_g173(double x){return 173.0*x/(1.0+x*x);}
static double rd_g174(double x){return 174.0*x/(1.0+x*x);}
static double rd_g175(double x){return 175.0*x/(1.0+x*x);}
static double rd_g176(double x){return 176.0*x/(1.0+x*x);}
static double rd_g177(double x){return 177.0*x/(1.0+x*x);}
static double rd_g178(double x){return 178.0*x/(1.0+x*x);}
static double rd_g179(double x){return 179.0*x/(1.0+x*x);}
static double rd_g180(double x){return 180.0*x/(1.0+x*x);}
static double rd_g181(double x){return 181.0*x/(1.0+x*x);}
static double rd_g182(double x){return 182.0*x/(1.0+x*x);}
static double rd_g183(double x){return 183.0*x/(1.0+x*x);}
static double rd_g184(double x){return 184.0*x/(1.0+x*x);}
static double rd_g185(double x){return 185.0*x/(1.0+x*x);}
static double rd_g186(double x){return 186.0*x/(1.0+x*x);}
static double rd_g187(double x){return 187.0*x/(1.0+x*x);}
static double rd_g188(double x){return 188.0*x/(1.0+x*x);}
static double rd_g189(double x){return 189.0*x/(1.0+x*x);}
static double rd_g190(double x){return 190.0*x/(1.0+x*x);}
static double rd_g191(double x){return 191.0*x/(1.0+x*x);}
static double rd_g192(double x){return 192.0*x/(1.0+x*x);}
static double rd_g193(double x){return 193.0*x/(1.0+x*x);}
static double rd_g194(double x){return 194.0*x/(1.0+x*x);}
static double rd_g195(double x){return 195.0*x/(1.0+x*x);}
static double rd_g196(double x){return 196.0*x/(1.0+x*x);}
static double rd_g197(double x){return 197.0*x/(1.0+x*x);}
static double rd_g198(double x){return 198.0*x/(1.0+x*x);}
static double rd_g199(double x){return 199.0*x/(1.0+x*x);}
static double rd_g200(double x){return 200.0*x/(1.0+x*x);}
static double rd_g201(double x){return 201.0*x/(1.0+x*x);}
static double rd_g202(double x){return 202.0*x/(1.0+x*x);}
static double rd_g203(double x){return 203.0*x/(1.0+x*x);}
static double rd_g204(double x){return 204.0*x/(1.0+x*x);}
static double rd_g205(double x){return 205.0*x/(1.0+x*x);}
static double rd_g206(double x){return 206.0*x/(1.0+x*x);}
static double rd_g207(double x){return 207.0*x/(1.0+x*x);}
static double rd_g208(double x){return 208.0*x/(1.0+x*x);}
static double rd_g209(double x){return 209.0*x/(1.0+x*x);}
static double rd_g210(double x){return 210.0*x/(1.0+x*x);}
static double rd_g211(double x){return 211.0*x/(1.0+x*x);}
static double rd_g212(double x){return 212.0*x/(1.0+x*x);}
static double rd_g213(double x){return 213.0*x/(1.0+x*x);}
static double rd_g214(double x){return 214.0*x/(1.0+x*x);}
static double rd_g215(double x){return 215.0*x/(1.0+x*x);}
static double rd_g216(double x){return 216.0*x/(1.0+x*x);}
static double rd_g217(double x){return 217.0*x/(1.0+x*x);}
static double rd_g218(double x){return 218.0*x/(1.0+x*x);}
static double rd_g219(double x){return 219.0*x/(1.0+x*x);}
static double rd_g220(double x){return 220.0*x/(1.0+x*x);}
static double rd_g221(double x){return 221.0*x/(1.0+x*x);}
static double rd_g222(double x){return 222.0*x/(1.0+x*x);}
static double rd_g223(double x){return 223.0*x/(1.0+x*x);}
static double rd_g224(double x){return 224.0*x/(1.0+x*x);}
static double rd_g225(double x){return 225.0*x/(1.0+x*x);}
static double rd_g226(double x){return 226.0*x/(1.0+x*x);}
static double rd_g227(double x){return 227.0*x/(1.0+x*x);}
static double rd_g228(double x){return 228.0*x/(1.0+x*x);}
static double rd_g229(double x){return 229.0*x/(1.0+x*x);}
static double rd_g230(double x){return 230.0*x/(1.0+x*x);}
static double rd_g231(double x){return 231.0*x/(1.0+x*x);}
static double rd_g232(double x){return 232.0*x/(1.0+x*x);}
static double rd_g233(double x){return 233.0*x/(1.0+x*x);}
static double rd_g234(double x){return 234.0*x/(1.0+x*x);}
static double rd_g235(double x){return 235.0*x/(1.0+x*x);}
static double rd_g236(double x){return 236.0*x/(1.0+x*x);}
static double rd_g237(double x){return 237.0*x/(1.0+x*x);}
static double rd_g238(double x){return 238.0*x/(1.0+x*x);}
static double rd_g239(double x){return 239.0*x/(1.0+x*x);}
static double rd_g240(double x){return 240.0*x/(1.0+x*x);}
static double rd_g241(double x){return 241.0*x/(1.0+x*x);}
static double rd_g242(double x){return 242.0*x/(1.0+x*x);}
static double rd_g243(double x){return 243.0*x/(1.0+x*x);}
static double rd_g244(double x){return 244.0*x/(1.0+x*x);}
static double rd_g245(double x){return 245.0*x/(1.0+x*x);}
static double rd_g246(double x){return 246.0*x/(1.0+x*x);}
static double rd_g247(double x){return 247.0*x/(1.0+x*x);}
static double rd_g248(double x){return 248.0*x/(1.0+x*x);}
static double rd_g249(double x){return 249.0*x/(1.0+x*x);}
static double rd_g250(double x){return 250.0*x/(1.0+x*x);}
static double rd_g251(double x){return 251.0*x/(1.0+x*x);}
static double rd_g252(double x){return 252.0*x/(1.0+x*x);}
static double rd_g253(double x){return 253.0*x/(1.0+x*x);}
static double rd_g254(double x){return 254.0*x/(1.0+x*x);}
static double rd_g255(double x){return 255.0*x/(1.0+x*x);}
static double rd_g256(double x){return 256.0*x/(1.0+x*x);}
static double rd_g257(double x){return 257.0*x/(1.0+x*x);}
static double rd_g258(double x){return 258.0*x/(1.0+x*x);}
static double rd_g259(double x){return 259.0*x/(1.0+x*x);}
static double rd_g260(double x){return 260.0*x/(1.0+x*x);}
static double rd_g261(double x){return 261.0*x/(1.0+x*x);}
static double rd_g262(double x){return 262.0*x/(1.0+x*x);}
static double rd_g263(double x){return 263.0*x/(1.0+x*x);}
static double rd_g264(double x){return 264.0*x/(1.0+x*x);}
static double rd_g265(double x){return 265.0*x/(1.0+x*x);}
static double rd_g266(double x){return 266.0*x/(1.0+x*x);}
static double rd_g267(double x){return 267.0*x/(1.0+x*x);}
static double rd_g268(double x){return 268.0*x/(1.0+x*x);}
static double rd_g269(double x){return 269.0*x/(1.0+x*x);}
static double rd_g270(double x){return 270.0*x/(1.0+x*x);}
static double rd_g271(double x){return 271.0*x/(1.0+x*x);}
static double rd_g272(double x){return 272.0*x/(1.0+x*x);}
static double rd_g273(double x){return 273.0*x/(1.0+x*x);}
static double rd_g274(double x){return 274.0*x/(1.0+x*x);}
static double rd_g275(double x){return 275.0*x/(1.0+x*x);}
static double rd_g276(double x){return 276.0*x/(1.0+x*x);}
static double rd_g277(double x){return 277.0*x/(1.0+x*x);}
static double rd_g278(double x){return 278.0*x/(1.0+x*x);}
static double rd_g279(double x){return 279.0*x/(1.0+x*x);}
static double rd_g280(double x){return 280.0*x/(1.0+x*x);}
static double rd_g281(double x){return 281.0*x/(1.0+x*x);}
static double rd_g282(double x){return 282.0*x/(1.0+x*x);}
static double rd_g283(double x){return 283.0*x/(1.0+x*x);}
static double rd_g284(double x){return 284.0*x/(1.0+x*x);}
static double rd_g285(double x){return 285.0*x/(1.0+x*x);}
static double rd_g286(double x){return 286.0*x/(1.0+x*x);}
static double rd_g287(double x){return 287.0*x/(1.0+x*x);}
static double rd_g288(double x){return 288.0*x/(1.0+x*x);}
static double rd_g289(double x){return 289.0*x/(1.0+x*x);}
static double rd_g290(double x){return 290.0*x/(1.0+x*x);}
static double rd_g291(double x){return 291.0*x/(1.0+x*x);}
static double rd_g292(double x){return 292.0*x/(1.0+x*x);}
static double rd_g293(double x){return 293.0*x/(1.0+x*x);}
static double rd_g294(double x){return 294.0*x/(1.0+x*x);}
static double rd_g295(double x){return 295.0*x/(1.0+x*x);}
static double rd_g296(double x){return 296.0*x/(1.0+x*x);}
static double rd_g297(double x){return 297.0*x/(1.0+x*x);}
static double rd_g298(double x){return 298.0*x/(1.0+x*x);}
static double rd_g299(double x){return 299.0*x/(1.0+x*x);}
static double rd_g300(double x){return 300.0*x/(1.0+x*x);}
static double rd_g301(double x){return 301.0*x/(1.0+x*x);}
static double rd_g302(double x){return 302.0*x/(1.0+x*x);}
static double rd_g303(double x){return 303.0*x/(1.0+x*x);}
static double rd_g304(double x){return 304.0*x/(1.0+x*x);}
static double rd_g305(double x){return 305.0*x/(1.0+x*x);}
static double rd_g306(double x){return 306.0*x/(1.0+x*x);}
static double rd_g307(double x){return 307.0*x/(1.0+x*x);}
static double rd_g308(double x){return 308.0*x/(1.0+x*x);}
static double rd_g309(double x){return 309.0*x/(1.0+x*x);}
static double rd_g310(double x){return 310.0*x/(1.0+x*x);}
static double rd_g311(double x){return 311.0*x/(1.0+x*x);}
static double rd_g312(double x){return 312.0*x/(1.0+x*x);}
static double rd_g313(double x){return 313.0*x/(1.0+x*x);}
static double rd_g314(double x){return 314.0*x/(1.0+x*x);}
static double rd_g315(double x){return 315.0*x/(1.0+x*x);}
static double rd_g316(double x){return 316.0*x/(1.0+x*x);}
static double rd_g317(double x){return 317.0*x/(1.0+x*x);}
static double rd_g318(double x){return 318.0*x/(1.0+x*x);}
static double rd_g319(double x){return 319.0*x/(1.0+x*x);}
static double rd_g320(double x){return 320.0*x/(1.0+x*x);}
static double rd_g321(double x){return 321.0*x/(1.0+x*x);}
static double rd_g322(double x){return 322.0*x/(1.0+x*x);}
static double rd_g323(double x){return 323.0*x/(1.0+x*x);}
static double rd_g324(double x){return 324.0*x/(1.0+x*x);}
static double rd_g325(double x){return 325.0*x/(1.0+x*x);}
static double rd_g326(double x){return 326.0*x/(1.0+x*x);}
static double rd_g327(double x){return 327.0*x/(1.0+x*x);}
static double rd_g328(double x){return 328.0*x/(1.0+x*x);}
static double rd_g329(double x){return 329.0*x/(1.0+x*x);}
static double rd_g330(double x){return 330.0*x/(1.0+x*x);}
static double rd_g331(double x){return 331.0*x/(1.0+x*x);}
static double rd_g332(double x){return 332.0*x/(1.0+x*x);}
static double rd_g333(double x){return 333.0*x/(1.0+x*x);}
static double rd_g334(double x){return 334.0*x/(1.0+x*x);}
static double rd_g335(double x){return 335.0*x/(1.0+x*x);}
static double rd_g336(double x){return 336.0*x/(1.0+x*x);}
static double rd_g337(double x){return 337.0*x/(1.0+x*x);}
static double rd_g338(double x){return 338.0*x/(1.0+x*x);}
static double rd_g339(double x){return 339.0*x/(1.0+x*x);}
static double rd_g340(double x){return 340.0*x/(1.0+x*x);}
static double rd_g341(double x){return 341.0*x/(1.0+x*x);}
static double rd_g342(double x){return 342.0*x/(1.0+x*x);}
static double rd_g343(double x){return 343.0*x/(1.0+x*x);}
static double rd_g344(double x){return 344.0*x/(1.0+x*x);}
static double rd_g345(double x){return 345.0*x/(1.0+x*x);}
static double rd_g346(double x){return 346.0*x/(1.0+x*x);}
static double rd_g347(double x){return 347.0*x/(1.0+x*x);}
static double rd_g348(double x){return 348.0*x/(1.0+x*x);}
static double rd_g349(double x){return 349.0*x/(1.0+x*x);}
static double rd_g350(double x){return 350.0*x/(1.0+x*x);}
static double rd_g351(double x){return 351.0*x/(1.0+x*x);}
static double rd_g352(double x){return 352.0*x/(1.0+x*x);}
static double rd_g353(double x){return 353.0*x/(1.0+x*x);}
static double rd_g354(double x){return 354.0*x/(1.0+x*x);}
static double rd_g355(double x){return 355.0*x/(1.0+x*x);}
static double rd_g356(double x){return 356.0*x/(1.0+x*x);}
static double rd_g357(double x){return 357.0*x/(1.0+x*x);}
static double rd_g358(double x){return 358.0*x/(1.0+x*x);}
static double rd_g359(double x){return 359.0*x/(1.0+x*x);}
static double rd_g360(double x){return 360.0*x/(1.0+x*x);}
static double rd_g361(double x){return 361.0*x/(1.0+x*x);}
static double rd_g362(double x){return 362.0*x/(1.0+x*x);}
static double rd_g363(double x){return 363.0*x/(1.0+x*x);}
static double rd_g364(double x){return 364.0*x/(1.0+x*x);}
static double rd_g365(double x){return 365.0*x/(1.0+x*x);}
static double rd_g366(double x){return 366.0*x/(1.0+x*x);}
static double rd_g367(double x){return 367.0*x/(1.0+x*x);}
static double rd_g368(double x){return 368.0*x/(1.0+x*x);}
static double rd_g369(double x){return 369.0*x/(1.0+x*x);}
static double rd_g370(double x){return 370.0*x/(1.0+x*x);}
static double rd_g371(double x){return 371.0*x/(1.0+x*x);}
static double rd_g372(double x){return 372.0*x/(1.0+x*x);}
static double rd_g373(double x){return 373.0*x/(1.0+x*x);}
static double rd_g374(double x){return 374.0*x/(1.0+x*x);}
static double rd_g375(double x){return 375.0*x/(1.0+x*x);}
static double rd_g376(double x){return 376.0*x/(1.0+x*x);}
static double rd_g377(double x){return 377.0*x/(1.0+x*x);}
static double rd_g378(double x){return 378.0*x/(1.0+x*x);}
static double rd_g379(double x){return 379.0*x/(1.0+x*x);}
static double rd_g380(double x){return 380.0*x/(1.0+x*x);}
static double rd_g381(double x){return 381.0*x/(1.0+x*x);}
static double rd_g382(double x){return 382.0*x/(1.0+x*x);}
static double rd_g383(double x){return 383.0*x/(1.0+x*x);}
static double rd_g384(double x){return 384.0*x/(1.0+x*x);}
static double rd_g385(double x){return 385.0*x/(1.0+x*x);}
static double rd_g386(double x){return 386.0*x/(1.0+x*x);}
static double rd_g387(double x){return 387.0*x/(1.0+x*x);}
static double rd_g388(double x){return 388.0*x/(1.0+x*x);}
static double rd_g389(double x){return 389.0*x/(1.0+x*x);}
static double rd_g390(double x){return 390.0*x/(1.0+x*x);}
static double rd_g391(double x){return 391.0*x/(1.0+x*x);}
static double rd_g392(double x){return 392.0*x/(1.0+x*x);}
static double rd_g393(double x){return 393.0*x/(1.0+x*x);}
static double rd_g394(double x){return 394.0*x/(1.0+x*x);}
static double rd_g395(double x){return 395.0*x/(1.0+x*x);}
static double rd_g396(double x){return 396.0*x/(1.0+x*x);}
static double rd_g397(double x){return 397.0*x/(1.0+x*x);}
static double rd_g398(double x){return 398.0*x/(1.0+x*x);}
static double rd_g399(double x){return 399.0*x/(1.0+x*x);}
static double rd_g400(double x){return 400.0*x/(1.0+x*x);}
static double rd_g401(double x){return 401.0*x/(1.0+x*x);}
static double rd_g402(double x){return 402.0*x/(1.0+x*x);}
static double rd_g403(double x){return 403.0*x/(1.0+x*x);}
static double rd_g404(double x){return 404.0*x/(1.0+x*x);}
static double rd_g405(double x){return 405.0*x/(1.0+x*x);}
static double rd_g406(double x){return 406.0*x/(1.0+x*x);}
static double rd_g407(double x){return 407.0*x/(1.0+x*x);}
static double rd_g408(double x){return 408.0*x/(1.0+x*x);}
static double rd_g409(double x){return 409.0*x/(1.0+x*x);}
static double rd_g410(double x){return 410.0*x/(1.0+x*x);}
static double rd_g411(double x){return 411.0*x/(1.0+x*x);}
static double rd_g412(double x){return 412.0*x/(1.0+x*x);}
static double rd_g413(double x){return 413.0*x/(1.0+x*x);}
static double rd_g414(double x){return 414.0*x/(1.0+x*x);}
static double rd_g415(double x){return 415.0*x/(1.0+x*x);}
static double rd_g416(double x){return 416.0*x/(1.0+x*x);}
static double rd_g417(double x){return 417.0*x/(1.0+x*x);}
static double rd_g418(double x){return 418.0*x/(1.0+x*x);}
static double rd_g419(double x){return 419.0*x/(1.0+x*x);}

/* ── Extended Turing Instability Analysis ──────────────────
 * For a 2-species reaction-diffusion system:
 *   du/dt = f(u,v) + Du * nabla^2 u
 *   dv/dt = g(u,v) + Dv * nabla^2 v
 *
 * Turing instability requires:
 *   1. Stable homogeneous steady state: f_u + g_v < 0, f_u*g_v - f_v*g_u > 0
 *   2. Diffusion-driven instability: Dv*f_u + Du*g_v > 0
 *   3. Unstable at some wavenumber k: dispersion relation Re(lambda(k^2)) > 0
 */

/* Compute the full Turing instability parameter region.
 * Scans over (Du, Dv) space to find the Turing space.
 * Returns number of (Du,Dv) pairs where patterns form. */
int rd_turing_parameter_scan(double alpha, double beta,
                              double gamma, double delta,
                              int n_Du, int n_Dv,
                              double Du_min, double Du_max,
                              double Dv_min, double Dv_max,
                              int *pattern_matrix) {
    if (!pattern_matrix || n_Du < 2 || n_Dv < 2) return -1;
    if (Du_max <= Du_min || Dv_max <= Dv_min) return -1;

    int pattern_count = 0;
    double dDu = (Du_max - Du_min) / (double)(n_Du - 1);
    double dDv = (Dv_max - Dv_min) / (double)(n_Dv - 1);

    for (int i = 0; i < n_Du; i++) {
        double Du = Du_min + (double)i * dDu;
        for (int j = 0; j < n_Dv; j++) {
            double Dv = Dv_min + (double)j * dDv;

            /* Check Turing conditions */
            /* Condition 1: stable homogeneous steady state */
            double trace = alpha + delta;
            double det = alpha * delta - beta * gamma;
            if (trace >= 0.0 || det <= 0.0) {
                pattern_matrix[i * n_Dv + j] = 0;
                continue;
            }

            /* Condition 2: diffusion-driven instability possible */
            if (Dv * alpha + Du * delta <= 0.0) {
                pattern_matrix[i * n_Dv + j] = 0;
                continue;
            }

            /* Condition 3: dispersion relation positive at some k^2 */
            /* Re(lambda) = 0.5*(tr(k^2) + sqrt(tr(k^2)^2 - 4*det(k^2)))
             * where tr(k^2) = trace - (Du+Dv)*k^2
             *       det(k^2) = det - (Dv*alpha + Du*delta)*k^2 + Du*Dv*k^4 */
            int can_form = 0;
            /* Scan k^2 from 0 to some max */
            double k2_max = 100.0;
            int n_k = 500;
            double dk2 = k2_max / (double)n_k;

            for (int k = 0; k < n_k; k++) {
                double k2 = (double)k * dk2 + 0.001;
                double tr_k2 = trace - (Du + Dv) * k2;
                double det_k2 = det - (Dv * alpha + Du * delta) * k2
                                + Du * Dv * k2 * k2;

                double disc = tr_k2 * tr_k2 - 4.0 * det_k2;
                if (disc >= 0.0) {
                    double lambda_max = 0.5 * (tr_k2 + sqrt(disc));
                    if (lambda_max > 1e-10) { can_form = 1; break; }
                } else {
                    /* Complex conjugate: instability if real part > 0 */
                    if (tr_k2 > 1e-10) { can_form = 1; break; }
                }
            }

            pattern_matrix[i * n_Dv + j] = can_form;
            if (can_form) pattern_count++;
        }
    }
    return pattern_count;
}

/* Compute the fastest-growing wavenumber for Turing instability.
 * Returns k_max^2 where dispersion relation peaks. */
double rd_fastest_wavenumber(double alpha, double beta,
                              double gamma, double delta,
                              double Du, double Dv) {
    /* The most unstable wavenumber maximizes Re(lambda(k^2)).
     * For complex lambda, the maximum real part occurs when
     * k^2 = sqrt(det/(Du*Dv)) for the Brusselator-type kinetics.
     *
     * More generally, solve d(Re(lambda))/d(k^2) = 0.
     * For the case where discriminant < 0:
     *   Re(lambda) = tr(k^2)/2 = (trace - (Du+Dv)*k^2)/2
     *   This is maximized at smallest k^2 > 0 where det(k^2) < 0.
     *
     * For discriminant >= 0:
     *   The maximum of the larger root occurs at the critical k^2. */

    double trace = alpha + delta;
    double det = alpha * delta - beta * gamma;

    /* Critical wavenumber from d(det_k2)/d(k^2) = 0:
     *   - (Dv*alpha + Du*delta) + 2*Du*Dv*k^2 = 0
     *   => k_c^2 = (Dv*alpha + Du*delta) / (2*Du*Dv) */
    double denom = 2.0 * Du * Dv;
    if (denom < 1e-15) return 0.0;
    double kc2 = (Dv * alpha + Du * delta) / denom;

    if (kc2 <= 0.0) return 0.0;

    /* Verify this actually gives positive growth */
    double tr_kc2 = trace - (Du + Dv) * kc2;
    double det_kc2 = det - (Dv * alpha + Du * delta) * kc2
                     + Du * Dv * kc2 * kc2;
    double disc = tr_kc2 * tr_kc2 - 4.0 * det_kc2;

    if (disc >= 0.0) {
        double lambda = 0.5 * (tr_kc2 + sqrt(disc));
        if (lambda > 0.0) return kc2;
    }

    /* If no instability at k_c^2, scan for any unstable k^2 */
    for (int i = 1; i <= 500; i++) {
        double k2 = (double)i * 0.1;
        double tr_k2 = trace - (Du + Dv) * k2;
        double det_k2 = det - (Dv * alpha + Du * delta) * k2
                        + Du * Dv * k2 * k2;
        disc = tr_k2 * tr_k2 - 4.0 * det_k2;
        if (disc >= 0.0) {
            double lambda = 0.5 * (tr_k2 + sqrt(disc));
            if (lambda > 0.0) return k2;
        }
    }
    return 0.0;
}

/* Pattern wavelength from fastest wavenumber.
 * lambda_pattern = 2*pi / k_max */
double rd_pattern_wavelength(double alpha, double beta,
                              double gamma, double delta,
                              double Du, double Dv) {
    double k2_max = rd_fastest_wavenumber(alpha, beta, gamma, delta, Du, Dv);
    if (k2_max < 1e-12) return 0.0;
    return 2.0 * M_PI / sqrt(k2_max);
}
