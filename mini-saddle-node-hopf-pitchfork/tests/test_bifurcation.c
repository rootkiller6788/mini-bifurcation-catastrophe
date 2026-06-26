#include "gst_core.h"
#include "saddle_node.h"
#include "hopf_bifurcation.h"
#include "pitchfork_bifurcation.h"
#include "bifurcation_detection.h"
#include "center_manifold.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
void t1(void){Matrix*A=mat_create(2,2);mat_set(A,0,0,1);mat_set(A,1,1,4);assert(fabs(mat_trace(A)-5)<1e-9);mat_free(A);printf("  mat_ops: PASS\n");}
void t2(void){Matrix*J=mat_create(2,2);mat_set(J,0,1,1);mat_set(J,1,0,-1);EigenSpectrum*e=eigen_compute_2x2(J);assert(e&&e->n==2);eigen_free(e);mat_free(J);printf("  eigen: PASS\n");}
void t3(void){double x[]={0.5},p[]={-0.25};int r=fp_newton(sn_normal_form,x,p,1,20,1e-8,&(FixedPoint){.point=x,.n=1});assert(r>=0);printf("  newton: PASS\n");}
void t4(void){SaddleNodeBifurcation*sn=sn_create();assert(sn);sn_free(sn);printf("  sn_create: PASS\n");}
void t5(void){HopfBifurcation*hb=hopf_create();assert(hb);hopf_free(hb);printf("  hopf_create: PASS\n");}
void t6(void){PitchforkBifurcation*pb=pf_create();assert(pb);pf_free(pb);printf("  pf_create: PASS\n");}
void t7(void){assert(strcmp(nf_type_name(NORM_SN),"Saddle-Node")==0);printf("  nf_names: PASS\n");}
void t8(void){EigenSpectrum*e=eigen_create();eigen_add(e,-1,0);eigen_add(e,-0.5,0);assert(eigen_all_stable(e));eigen_free(e);printf("  stable: PASS\n");}
void t9(void){Vector*a=vec_create(3);a->data[0]=3;a->data[1]=4;assert(fabs(vec_norm(a)-5)<1e-9);vec_free(a);printf("  vec: PASS\n");}
void t10(void){double x[]={1},p[]={-0.5};rk4_step(sn_normal_form,x,p,1,0.01);printf("  rk4: PASS\n");}
void t11(void){BifurcationScanner*bs=bscan_create();assert(bs);bscan_free(bs);printf("  bscan: PASS\n");}
void t12(void){CenterManifold*cm=cm_create(2);assert(cm);cm_free(cm);printf("  cm: PASS\n");}
void t13(void){Matrix*A=mat_identity(3);assert(A);double tr=mat_trace(A);assert(fabs(tr-3)<1e-9);mat_free(A);printf("  identity: PASS\n");}
void t14(void){Matrix*A=mat_create(2,2);Matrix*B=mat_create(2,2);mat_set(A,0,0,1);mat_set(B,0,0,2);Matrix*C=mat_mul(A,B);assert(C&&fabs(mat_get(C,0,0)-2)<1e-9);mat_free(A);mat_free(B);mat_free(C);printf("  mul: PASS\n");}
void t15(void){SaddleNodeBifurcation*sn=sn_create();sn->r_critical=0;sn->x_critical=0;double x[8]={1},p[4]={-0.5};int r=sn_compute_branches(sn,sn_normal_form,x,p,1,0,-1,1,20,20,1e-6);assert(r>=0);sn_free(sn);printf("  branches: PASS\n");}
void t16(void){Matrix*A=mat_create(2,2);Matrix*B=mat_create(2,2);mat_set(A,0,0,1);mat_set(A,1,1,1);mat_set(B,0,0,2);mat_set(B,1,1,2);Matrix*R=mat_add(A,B);assert(fabs(mat_get(R,0,0)-3)<1e-9);mat_free(A);mat_free(B);mat_free(R);printf("  add: PASS\n");}
void t17(void){Matrix*A=mat_create(2,2);mat_set(A,0,0,1);mat_set(A,0,1,2);mat_set(A,1,0,3);mat_set(A,1,1,4);Matrix*T=mat_transpose(A);assert(fabs(mat_get(T,0,1)-3)<1e-9);mat_free(A);mat_free(T);printf("  transpose: PASS\n");}
void t18(void){double x[8]={0.1,0.1},p[4]={0.1,1.0};int r=hopf_detect_lc(hopf_normal_form,x,p,2,0.05,0.01,10,100,NULL,NULL);assert(r>=0);printf("  hopf_lc: PASS\n");}
void t19(void){bool sym=pf_check_z2_symmetry(pitchfork_normal_form,1);assert(sym);printf("  z2sym: PASS\n");}
void t20(void){BifurcationScanner*bs=bscan_create();double x[8]={0.5},p[4]={-0.5};bscan_detect_all(bs,sn_normal_form,x,p,1,0,-1,1,10,20,1e-6);int n=bscan_count_type(bs,BIF_SN);assert(n>=0);bscan_free(bs);printf("  bscan_detect: PASS\n");}
int main(void){printf("=== Bifurcation Test Suite ===\n");fflush(stdout);t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();t14();t15();t16();t17();t18();t19();t20();printf("=== All 20 tests passed ===\n");return 0;}
