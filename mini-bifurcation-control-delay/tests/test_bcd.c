/* test_bcd.c */
#include "bcd_core.h"
#include "bcd_hopf.h"
#include "bcd_control.h"
#include "bcd_stability.h"
#include "bcd_washout.h"
#include "bcd_amplitude.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
static void td(int n,const double* x,const double* xd,double t,double* dx,void* p){(void)t;(void)p;dx[0]=-x[0]+0.5*xd[0];if(n>1)dx[1]=-2*x[1]+0.3*xd[1];}
int main(void){
assert(1+1==2); assert(BCD_PI>3.0); assert(BCD_E>2.0); printf("=== BCD Tests ===\n");
BCDDDE* d=bcd_dde_create(2,1.0,1000);assert(d);assert(d->dim==2);
bcd_dde_set_history(d,(double[]){0.1,0.2},2);assert(bcd_dde_step(d,td,NULL,0.01)==0);
assert(bcd_dde_get_delayed(d,0)>=0);
double out[100];bcd_dde_simulate(d,td,NULL,1.0,0.01,out,50);
BCDDFControl* c=bcd_control_create(0.5,1.0,2);assert(c);double u[2];assert(bcd_control_apply(c,(double[]){0.1,0.2},0,u)==0);
bcd_control_update(c,(double[]){0.1,0.2});bcd_control_free(c);
BCDWashoutFilter* w=bcd_washout_create(2.0,0.5,2);assert(w);double y[2];assert(bcd_washout_step(w,(double[]){0.1,0.2},0.01,y)==0);bcd_washout_free(w);
BCDHopfResult* hr=bcd_hopf_detect(d,td,NULL,0.5,2.0,20);assert(hr);bcd_hopf_result_free(hr);
BCDHopfBifPoint hbp=bcd_hopf_normal_form(d,td,NULL,1.0);assert(hbp.mu_crit>0);
BCDCoupledDDE* cd=bcd_coupled_create(3,2,0.5,1.0);assert(cd);bcd_coupled_free(cd);
double w0=bcd_lambert_w(1.0,0);assert(w0>0.4&&w0<0.6);
double A[]={-1,0,0,-2},B[]={0.5,0,0,0.3};BCDDelaySystem* ds=bcd_delay_system_create(A,B,2,1.0);assert(ds);
bcd_delay_system_eigenvalues(ds,10);assert(ds->n_eig==10);assert(bcd_num_unstable_roots(ds)>=0);
double sr=bcd_stability_radius(ds);assert(sr>0);bcd_delay_system_free(ds);
bcd_dde_free(NULL);bcd_control_free(NULL);bcd_washout_free(NULL);bcd_hopf_result_free(NULL);
assert(bcd_lambert_w(0.0,0) < 1e-10);
BCDHopfResult* hr2 = bcd_hopf_result_create(); assert(hr2); bcd_hopf_result_free(hr2);
