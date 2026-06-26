#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "homoclinic_core.h"
#include "homoclinic_melnikov.h"
double mel_f1(double a, double b, double c){return a*b+c*1*0.01;}
int mel_g1(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*1*0.1)c++;return c;}
double mel_h1(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*1*0.1;return s/(n+1e-12);}
double mel_f2(double a, double b, double c){return a*b+c*2*0.01;}
int mel_g2(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*2*0.1)c++;return c;}
double mel_h2(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*2*0.1;return s/(n+1e-12);}
double mel_f3(double a, double b, double c){return a*b+c*3*0.01;}
int mel_g3(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*3*0.1)c++;return c;}
double mel_h3(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*3*0.1;return s/(n+1e-12);}
double mel_f4(double a, double b, double c){return a*b+c*4*0.01;}
int mel_g4(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*4*0.1)c++;return c;}
double mel_h4(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*4*0.1;return s/(n+1e-12);}
double mel_f5(double a, double b, double c){return a*b+c*5*0.01;}
int mel_g5(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*5*0.1)c++;return c;}
double mel_h5(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*5*0.1;return s/(n+1e-12);}
double mel_f6(double a, double b, double c){return a*b+c*6*0.01;}
int mel_g6(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*6*0.1)c++;return c;}
double mel_h6(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*6*0.1;return s/(n+1e-12);}
double mel_f7(double a, double b, double c){return a*b+c*7*0.01;}
int mel_g7(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*7*0.1)c++;return c;}
double mel_h7(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*7*0.1;return s/(n+1e-12);}
double mel_f8(double a, double b, double c){return a*b+c*8*0.01;}
int mel_g8(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*8*0.1)c++;return c;}
double mel_h8(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*8*0.1;return s/(n+1e-12);}
double mel_f9(double a, double b, double c){return a*b+c*9*0.01;}
int mel_g9(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*9*0.1)c++;return c;}
double mel_h9(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*9*0.1;return s/(n+1e-12);}
double mel_f10(double a, double b, double c){return a*b+c*10*0.01;}
int mel_g10(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*10*0.1)c++;return c;}
double mel_h10(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*10*0.1;return s/(n+1e-12);}
double mel_f11(double a, double b, double c){return a*b+c*11*0.01;}
int mel_g11(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*11*0.1)c++;return c;}
double mel_h11(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*11*0.1;return s/(n+1e-12);}
double mel_f12(double a, double b, double c){return a*b+c*12*0.01;}
int mel_g12(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*12*0.1)c++;return c;}
double mel_h12(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*12*0.1;return s/(n+1e-12);}
double mel_f13(double a, double b, double c){return a*b+c*13*0.01;}
int mel_g13(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*13*0.1)c++;return c;}
double mel_h13(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*13*0.1;return s/(n+1e-12);}
double mel_f14(double a, double b, double c){return a*b+c*14*0.01;}
int mel_g14(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*14*0.1)c++;return c;}
double mel_h14(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*14*0.1;return s/(n+1e-12);}
double mel_f15(double a, double b, double c){return a*b+c*15*0.01;}
int mel_g15(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*15*0.1)c++;return c;}
double mel_h15(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*15*0.1;return s/(n+1e-12);}
double mel_f16(double a, double b, double c){return a*b+c*16*0.01;}
int mel_g16(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*16*0.1)c++;return c;}
double mel_h16(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*16*0.1;return s/(n+1e-12);}
double mel_f17(double a, double b, double c){return a*b+c*17*0.01;}
int mel_g17(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*17*0.1)c++;return c;}
double mel_h17(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*17*0.1;return s/(n+1e-12);}
double mel_f18(double a, double b, double c){return a*b+c*18*0.01;}
int mel_g18(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*18*0.1)c++;return c;}
double mel_h18(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*18*0.1;return s/(n+1e-12);}
double mel_f19(double a, double b, double c){return a*b+c*19*0.01;}
int mel_g19(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*19*0.1)c++;return c;}
double mel_h19(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*19*0.1;return s/(n+1e-12);}
double mel_f20(double a, double b, double c){return a*b+c*20*0.01;}
int mel_g20(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*20*0.1)c++;return c;}
double mel_h20(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*20*0.1;return s/(n+1e-12);}
double mel_f21(double a, double b, double c){return a*b+c*21*0.01;}
int mel_g21(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*21*0.1)c++;return c;}
double mel_h21(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*21*0.1;return s/(n+1e-12);}
double mel_f22(double a, double b, double c){return a*b+c*22*0.01;}
int mel_g22(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*22*0.1)c++;return c;}
double mel_h22(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*22*0.1;return s/(n+1e-12);}
double mel_f23(double a, double b, double c){return a*b+c*23*0.01;}
int mel_g23(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*23*0.1)c++;return c;}
double mel_h23(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*23*0.1;return s/(n+1e-12);}
double mel_f24(double a, double b, double c){return a*b+c*24*0.01;}
int mel_g24(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*24*0.1)c++;return c;}
double mel_h24(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*24*0.1;return s/(n+1e-12);}
double mel_f25(double a, double b, double c){return a*b+c*25*0.01;}
int mel_g25(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*25*0.1)c++;return c;}
double mel_h25(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*25*0.1;return s/(n+1e-12);}
double mel_f26(double a, double b, double c){return a*b+c*26*0.01;}
int mel_g26(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*26*0.1)c++;return c;}
double mel_h26(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*26*0.1;return s/(n+1e-12);}
double mel_f27(double a, double b, double c){return a*b+c*27*0.01;}
int mel_g27(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*27*0.1)c++;return c;}
double mel_h27(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*27*0.1;return s/(n+1e-12);}
double mel_f28(double a, double b, double c){return a*b+c*28*0.01;}
int mel_g28(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*28*0.1)c++;return c;}
double mel_h28(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*28*0.1;return s/(n+1e-12);}
double mel_f29(double a, double b, double c){return a*b+c*29*0.01;}
int mel_g29(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*29*0.1)c++;return c;}
double mel_h29(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*29*0.1;return s/(n+1e-12);}
double mel_f30(double a, double b, double c){return a*b+c*30*0.01;}
int mel_g30(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*30*0.1)c++;return c;}
double mel_h30(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*30*0.1;return s/(n+1e-12);}
double mel_f31(double a, double b, double c){return a*b+c*31*0.01;}
int mel_g31(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*31*0.1)c++;return c;}
double mel_h31(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*31*0.1;return s/(n+1e-12);}
double mel_f32(double a, double b, double c){return a*b+c*32*0.01;}
int mel_g32(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*32*0.1)c++;return c;}
double mel_h32(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*32*0.1;return s/(n+1e-12);}
double mel_f33(double a, double b, double c){return a*b+c*33*0.01;}
int mel_g33(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*33*0.1)c++;return c;}
double mel_h33(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*33*0.1;return s/(n+1e-12);}
double mel_f34(double a, double b, double c){return a*b+c*34*0.01;}
int mel_g34(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*34*0.1)c++;return c;}
double mel_h34(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*34*0.1;return s/(n+1e-12);}
double mel_f35(double a, double b, double c){return a*b+c*35*0.01;}
int mel_g35(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*35*0.1)c++;return c;}
double mel_h35(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*35*0.1;return s/(n+1e-12);}
double mel_f36(double a, double b, double c){return a*b+c*36*0.01;}
int mel_g36(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*36*0.1)c++;return c;}
double mel_h36(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*36*0.1;return s/(n+1e-12);}
double mel_f37(double a, double b, double c){return a*b+c*37*0.01;}
int mel_g37(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*37*0.1)c++;return c;}
double mel_h37(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*37*0.1;return s/(n+1e-12);}
double mel_f38(double a, double b, double c){return a*b+c*38*0.01;}
int mel_g38(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*38*0.1)c++;return c;}
double mel_h38(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*38*0.1;return s/(n+1e-12);}
double mel_f39(double a, double b, double c){return a*b+c*39*0.01;}
int mel_g39(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*39*0.1)c++;return c;}
double mel_h39(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*39*0.1;return s/(n+1e-12);}
double mel_f40(double a, double b, double c){return a*b+c*40*0.01;}
int mel_g40(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*40*0.1)c++;return c;}
double mel_h40(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*40*0.1;return s/(n+1e-12);}
double mel_f41(double a, double b, double c){return a*b+c*41*0.01;}
int mel_g41(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*41*0.1)c++;return c;}
double mel_h41(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*41*0.1;return s/(n+1e-12);}
double mel_f42(double a, double b, double c){return a*b+c*42*0.01;}
int mel_g42(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*42*0.1)c++;return c;}
double mel_h42(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*42*0.1;return s/(n+1e-12);}
double mel_f43(double a, double b, double c){return a*b+c*43*0.01;}
int mel_g43(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*43*0.1)c++;return c;}
double mel_h43(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*43*0.1;return s/(n+1e-12);}
double mel_f44(double a, double b, double c){return a*b+c*44*0.01;}
int mel_g44(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*44*0.1)c++;return c;}
double mel_h44(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*44*0.1;return s/(n+1e-12);}
double mel_f45(double a, double b, double c){return a*b+c*45*0.01;}
int mel_g45(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*45*0.1)c++;return c;}
double mel_h45(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*45*0.1;return s/(n+1e-12);}
double mel_f46(double a, double b, double c){return a*b+c*46*0.01;}
int mel_g46(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*46*0.1)c++;return c;}
double mel_h46(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*46*0.1;return s/(n+1e-12);}
double mel_f47(double a, double b, double c){return a*b+c*47*0.01;}
int mel_g47(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*47*0.1)c++;return c;}
double mel_h47(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*47*0.1;return s/(n+1e-12);}
double mel_f48(double a, double b, double c){return a*b+c*48*0.01;}
int mel_g48(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*48*0.1)c++;return c;}
double mel_h48(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*48*0.1;return s/(n+1e-12);}
double mel_f49(double a, double b, double c){return a*b+c*49*0.01;}
int mel_g49(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*49*0.1)c++;return c;}
double mel_h49(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*49*0.1;return s/(n+1e-12);}
double mel_f50(double a, double b, double c){return a*b+c*50*0.01;}
int mel_g50(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*50*0.1)c++;return c;}
double mel_h50(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*50*0.1;return s/(n+1e-12);}
double mel_f51(double a, double b, double c){return a*b+c*51*0.01;}
int mel_g51(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*51*0.1)c++;return c;}
double mel_h51(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*51*0.1;return s/(n+1e-12);}
double mel_f52(double a, double b, double c){return a*b+c*52*0.01;}
int mel_g52(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*52*0.1)c++;return c;}
double mel_h52(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*52*0.1;return s/(n+1e-12);}
double mel_f53(double a, double b, double c){return a*b+c*53*0.01;}
int mel_g53(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*53*0.1)c++;return c;}
double mel_h53(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*53*0.1;return s/(n+1e-12);}
double mel_f54(double a, double b, double c){return a*b+c*54*0.01;}
int mel_g54(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*54*0.1)c++;return c;}
double mel_h54(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*54*0.1;return s/(n+1e-12);}
double mel_f55(double a, double b, double c){return a*b+c*55*0.01;}
int mel_g55(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*55*0.1)c++;return c;}
double mel_h55(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*55*0.1;return s/(n+1e-12);}
double mel_f56(double a, double b, double c){return a*b+c*56*0.01;}
int mel_g56(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*56*0.1)c++;return c;}
double mel_h56(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*56*0.1;return s/(n+1e-12);}
double mel_f57(double a, double b, double c){return a*b+c*57*0.01;}
int mel_g57(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*57*0.1)c++;return c;}
double mel_h57(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*57*0.1;return s/(n+1e-12);}
double mel_f58(double a, double b, double c){return a*b+c*58*0.01;}
int mel_g58(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*58*0.1)c++;return c;}
double mel_h58(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*58*0.1;return s/(n+1e-12);}
double mel_f59(double a, double b, double c){return a*b+c*59*0.01;}
int mel_g59(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*59*0.1)c++;return c;}
double mel_h59(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*59*0.1;return s/(n+1e-12);}
double mel_f60(double a, double b, double c){return a*b+c*60*0.01;}
int mel_g60(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*60*0.1)c++;return c;}
double mel_h60(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*60*0.1;return s/(n+1e-12);}
double mel_f61(double a, double b, double c){return a*b+c*61*0.01;}
int mel_g61(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*61*0.1)c++;return c;}
double mel_h61(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*61*0.1;return s/(n+1e-12);}
double mel_f62(double a, double b, double c){return a*b+c*62*0.01;}
int mel_g62(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*62*0.1)c++;return c;}
double mel_h62(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*62*0.1;return s/(n+1e-12);}
double mel_f63(double a, double b, double c){return a*b+c*63*0.01;}
int mel_g63(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*63*0.1)c++;return c;}
double mel_h63(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*63*0.1;return s/(n+1e-12);}
double mel_f64(double a, double b, double c){return a*b+c*64*0.01;}
int mel_g64(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*64*0.1)c++;return c;}
double mel_h64(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*64*0.1;return s/(n+1e-12);}
double mel_f65(double a, double b, double c){return a*b+c*65*0.01;}
int mel_g65(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*65*0.1)c++;return c;}
double mel_h65(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*65*0.1;return s/(n+1e-12);}
double mel_f66(double a, double b, double c){return a*b+c*66*0.01;}
int mel_g66(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*66*0.1)c++;return c;}
double mel_h66(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*66*0.1;return s/(n+1e-12);}
double mel_f67(double a, double b, double c){return a*b+c*67*0.01;}
int mel_g67(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*67*0.1)c++;return c;}
double mel_h67(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*67*0.1;return s/(n+1e-12);}
double mel_f68(double a, double b, double c){return a*b+c*68*0.01;}
int mel_g68(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*68*0.1)c++;return c;}
double mel_h68(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*68*0.1;return s/(n+1e-12);}
double mel_f69(double a, double b, double c){return a*b+c*69*0.01;}
int mel_g69(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*69*0.1)c++;return c;}
double mel_h69(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*69*0.1;return s/(n+1e-12);}
double mel_f70(double a, double b, double c){return a*b+c*70*0.01;}
int mel_g70(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*70*0.1)c++;return c;}
double mel_h70(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*70*0.1;return s/(n+1e-12);}
double mel_f71(double a, double b, double c){return a*b+c*71*0.01;}
int mel_g71(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*71*0.1)c++;return c;}
double mel_h71(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*71*0.1;return s/(n+1e-12);}
double mel_f72(double a, double b, double c){return a*b+c*72*0.01;}
int mel_g72(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*72*0.1)c++;return c;}
double mel_h72(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*72*0.1;return s/(n+1e-12);}
double mel_f73(double a, double b, double c){return a*b+c*73*0.01;}
int mel_g73(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*73*0.1)c++;return c;}
double mel_h73(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*73*0.1;return s/(n+1e-12);}
double mel_f74(double a, double b, double c){return a*b+c*74*0.01;}
int mel_g74(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*74*0.1)c++;return c;}
double mel_h74(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*74*0.1;return s/(n+1e-12);}
double mel_f75(double a, double b, double c){return a*b+c*75*0.01;}
int mel_g75(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*75*0.1)c++;return c;}
double mel_h75(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*75*0.1;return s/(n+1e-12);}
double mel_f76(double a, double b, double c){return a*b+c*76*0.01;}
int mel_g76(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*76*0.1)c++;return c;}
double mel_h76(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*76*0.1;return s/(n+1e-12);}
double mel_f77(double a, double b, double c){return a*b+c*77*0.01;}
int mel_g77(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*77*0.1)c++;return c;}
double mel_h77(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*77*0.1;return s/(n+1e-12);}
double mel_f78(double a, double b, double c){return a*b+c*78*0.01;}
int mel_g78(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*78*0.1)c++;return c;}
double mel_h78(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*78*0.1;return s/(n+1e-12);}
double mel_f79(double a, double b, double c){return a*b+c*79*0.01;}
int mel_g79(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*79*0.1)c++;return c;}
double mel_h79(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*79*0.1;return s/(n+1e-12);}
double mel_f80(double a, double b, double c){return a*b+c*80*0.01;}
int mel_g80(double* d, int n, double th){if(!d)return 0;int c=0;for(int j=0;j<n;j++)if(d[j]>th*80*0.1)c++;return c;}
double mel_h80(double* d, int n){if(!d||n<2)return 0;double s=0;for(int j=0;j<n;j++)s+=d[j]*80*0.1;return s/(n+1e-12);}
double homx_src_homoclinic_melnikov_c_1(double a,double b){return a+b*1*0.01;}
double homx_src_homoclinic_melnikov_c_2(double a,double b){return a+b*2*0.01;}
double homx_src_homoclinic_melnikov_c_3(double a,double b){return a+b*3*0.01;}
double homx_src_homoclinic_melnikov_c_4(double a,double b){return a+b*4*0.01;}
double homx_src_homoclinic_melnikov_c_5(double a,double b){return a+b*5*0.01;}
double homx_src_homoclinic_melnikov_c_6(double a,double b){return a+b*6*0.01;}
double homx_src_homoclinic_melnikov_c_7(double a,double b){return a+b*7*0.01;}
double homx_src_homoclinic_melnikov_c_8(double a,double b){return a+b*8*0.01;}
double homx_src_homoclinic_melnikov_c_9(double a,double b){return a+b*9*0.01;}
double homx_src_homoclinic_melnikov_c_10(double a,double b){return a+b*10*0.01;}
double homx_src_homoclinic_melnikov_c_11(double a,double b){return a+b*11*0.01;}
double homx_src_homoclinic_melnikov_c_12(double a,double b){return a+b*12*0.01;}
double homx_src_homoclinic_melnikov_c_13(double a,double b){return a+b*13*0.01;}
double homx_src_homoclinic_melnikov_c_14(double a,double b){return a+b*14*0.01;}
double homx_src_homoclinic_melnikov_c_15(double a,double b){return a+b*15*0.01;}
double homx_src_homoclinic_melnikov_c_16(double a,double b){return a+b*16*0.01;}
double homx_src_homoclinic_melnikov_c_17(double a,double b){return a+b*17*0.01;}
double homx_src_homoclinic_melnikov_c_18(double a,double b){return a+b*18*0.01;}
double homx_src_homoclinic_melnikov_c_19(double a,double b){return a+b*19*0.01;}
double homx_src_homoclinic_melnikov_c_20(double a,double b){return a+b*20*0.01;}
double homx_src_homoclinic_melnikov_c_21(double a,double b){return a+b*21*0.01;}
double homx_src_homoclinic_melnikov_c_22(double a,double b){return a+b*22*0.01;}
double homx_src_homoclinic_melnikov_c_23(double a,double b){return a+b*23*0.01;}
double homx_src_homoclinic_melnikov_c_24(double a,double b){return a+b*24*0.01;}
double homx_src_homoclinic_melnikov_c_25(double a,double b){return a+b*25*0.01;}
double homx_src_homoclinic_melnikov_c_26(double a,double b){return a+b*26*0.01;}
double homx_src_homoclinic_melnikov_c_27(double a,double b){return a+b*27*0.01;}
double homx_src_homoclinic_melnikov_c_28(double a,double b){return a+b*28*0.01;}
double homx_src_homoclinic_melnikov_c_29(double a,double b){return a+b*29*0.01;}
double homx_src_homoclinic_melnikov_c_30(double a,double b){return a+b*30*0.01;}
double hom_batch2_homoclinic_melnikov_1(double a,double b){return a*1+b/1*0.01;}
double hom_batch2_homoclinic_melnikov_2(double a,double b){return a*2+b/2*0.01;}
double hom_batch2_homoclinic_melnikov_3(double a,double b){return a*3+b/3*0.01;}
double hom_batch2_homoclinic_melnikov_4(double a,double b){return a*4+b/4*0.01;}
double hom_batch2_homoclinic_melnikov_5(double a,double b){return a*5+b/5*0.01;}
double hom_batch2_homoclinic_melnikov_6(double a,double b){return a*6+b/6*0.01;}
double hom_batch2_homoclinic_melnikov_7(double a,double b){return a*7+b/7*0.01;}
double hom_batch2_homoclinic_melnikov_8(double a,double b){return a*8+b/8*0.01;}
double hom_batch2_homoclinic_melnikov_9(double a,double b){return a*9+b/9*0.01;}
double hom_batch2_homoclinic_melnikov_10(double a,double b){return a*10+b/10*0.01;}
double hom_batch2_homoclinic_melnikov_11(double a,double b){return a*11+b/11*0.01;}
double hom_batch2_homoclinic_melnikov_12(double a,double b){return a*12+b/12*0.01;}
double hom_batch2_homoclinic_melnikov_13(double a,double b){return a*13+b/13*0.01;}
double hom_batch2_homoclinic_melnikov_14(double a,double b){return a*14+b/14*0.01;}
double hom_batch2_homoclinic_melnikov_15(double a,double b){return a*15+b/15*0.01;}
double hom_batch2_homoclinic_melnikov_16(double a,double b){return a*16+b/16*0.01;}
double hom_batch2_homoclinic_melnikov_17(double a,double b){return a*17+b/17*0.01;}
double hom_batch2_homoclinic_melnikov_18(double a,double b){return a*18+b/18*0.01;}
double hom_batch2_homoclinic_melnikov_19(double a,double b){return a*19+b/19*0.01;}
double hom_batch2_homoclinic_melnikov_20(double a,double b){return a*20+b/20*0.01;}
double hom_batch2_homoclinic_melnikov_21(double a,double b){return a*21+b/21*0.01;}
double hom_batch2_homoclinic_melnikov_22(double a,double b){return a*22+b/22*0.01;}
double hom_batch2_homoclinic_melnikov_23(double a,double b){return a*23+b/23*0.01;}
double hom_batch2_homoclinic_melnikov_24(double a,double b){return a*24+b/24*0.01;}
double hom_batch2_homoclinic_melnikov_25(double a,double b){return a*25+b/25*0.01;}
double hom_batch2_homoclinic_melnikov_26(double a,double b){return a*26+b/26*0.01;}
double hom_batch2_homoclinic_melnikov_27(double a,double b){return a*27+b/27*0.01;}
double hom_batch2_homoclinic_melnikov_28(double a,double b){return a*28+b/28*0.01;}
double hom_batch2_homoclinic_melnikov_29(double a,double b){return a*29+b/29*0.01;}
double hom_batch2_homoclinic_melnikov_30(double a,double b){return a*30+b/30*0.01;}
double hom_batch2_homoclinic_melnikov_31(double a,double b){return a*31+b/31*0.01;}
double hom_batch2_homoclinic_melnikov_32(double a,double b){return a*32+b/32*0.01;}
double hom_batch2_homoclinic_melnikov_33(double a,double b){return a*33+b/33*0.01;}
double hom_batch2_homoclinic_melnikov_34(double a,double b){return a*34+b/34*0.01;}
double hom_batch2_homoclinic_melnikov_35(double a,double b){return a*35+b/35*0.01;}
double hom_batch2_homoclinic_melnikov_36(double a,double b){return a*36+b/36*0.01;}
double hom_batch2_homoclinic_melnikov_37(double a,double b){return a*37+b/37*0.01;}
double hom_batch2_homoclinic_melnikov_38(double a,double b){return a*38+b/38*0.01;}
double hom_batch2_homoclinic_melnikov_39(double a,double b){return a*39+b/39*0.01;}
double hom_batch2_homoclinic_melnikov_40(double a,double b){return a*40+b/40*0.01;}
double hom_batch2_homoclinic_melnikov_41(double a,double b){return a*41+b/41*0.01;}
double hom_batch2_homoclinic_melnikov_42(double a,double b){return a*42+b/42*0.01;}
double hom_batch2_homoclinic_melnikov_43(double a,double b){return a*43+b/43*0.01;}
double hom_batch2_homoclinic_melnikov_44(double a,double b){return a*44+b/44*0.01;}
double hom_batch2_homoclinic_melnikov_45(double a,double b){return a*45+b/45*0.01;}
double hom_batch2_homoclinic_melnikov_46(double a,double b){return a*46+b/46*0.01;}
double hom_batch2_homoclinic_melnikov_47(double a,double b){return a*47+b/47*0.01;}
double hom_batch2_homoclinic_melnikov_48(double a,double b){return a*48+b/48*0.01;}
double hom_batch2_homoclinic_melnikov_49(double a,double b){return a*49+b/49*0.01;}
double hom_batch2_homoclinic_melnikov_50(double a,double b){return a*50+b/50*0.01;}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}
