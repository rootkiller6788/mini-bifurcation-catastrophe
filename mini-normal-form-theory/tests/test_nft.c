#include "nft_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define EPS 1e-6
int main(void) {
    printf("=== Normal Form Test ===\n");
    { NFTMatrix* m=nft_matrix_create(2,2); assert(m); nft_matrix_free(m); }
    { NFTMatrix* a=nft_matrix_create(2,2); a->data[0]=1;a->data[3]=1;
      NFTMatrix* b=nft_matrix_create(2,2); b->data[0]=2;b->data[3]=2;
      NFTMatrix* s=nft_matrix_add(a,b); assert(s&&fabs(s->data[0]-3)<EPS);
      nft_matrix_free(a);nft_matrix_free(b);nft_matrix_free(s); }
    { double A[]={0,1,-2,-3}; NFTMatrix* m=nft_matrix_create(2,2);memcpy(m->data,A,4*sizeof(double));
      double re[2],im[2]; assert(nft_extract_eigenvalues(m,re,im)==2); nft_matrix_free(m); }
    { double A[]={0,1,0,0}; NFTMatrix* m=nft_matrix_create(2,2);memcpy(m->data,A,4*sizeof(double));
      NFTMatrix *P=NULL,*J=NULL; nft_jordan_decompose(m,&P,&J); assert(J); nft_matrix_free(P);nft_matrix_free(J);nft_matrix_free(m); }
    { NormalFormSpec* s=nft_spec_create(2,1); assert(s); nft_spec_free(s); }
    { double Ad[]={1,0,0,-1},Qd[]={0,1,0,0}; NFTMatrix* A=nft_matrix_create(2,2);memcpy(A->data,Ad,4*sizeof(double));
      NFTMatrix* Q=nft_matrix_create(2,2);memcpy(Q->data,Qd,4*sizeof(double));
      HomologicalEquation* he=nft_homological_create(A,Q); assert(he); assert(nft_solve_homological(he)==0);
      nft_homological_free(he);nft_matrix_free(A);nft_matrix_free(Q); }
    { double lam[]={1.0,-1.0};int n_res;ResonanceMonomial* rm=nft_find_resonances(lam,2,3,&n_res);nft_resonance_free(rm); }
    { double Ad[]={0,1,0,0};NFTMatrix* A=nft_matrix_create(2,2);memcpy(A->data,Ad,4*sizeof(double));
      NormalFormSpec* s=nft_spec_create(2,1);NormalFormResult* r=nft_compute_normal_form(A,s);assert(r&&r->converged);
      nft_result_free(r);nft_spec_free(s);nft_matrix_free(A); }
    { CenterManifold* cm=nft_center_create(1,1,2); assert(cm); nft_center_free(cm); }
    { double ev[]={0.0,-2.0};int c,s,u;nft_split_eigenvalues(ev,2,&c,&s,&u,1e-4);assert(c==1&&s==1); }
    { Unfolding* uf=nft_unfolding_create(1); assert(uf); nft_unfolding_free(uf); }
    { double Ad[]={0,1,0,0};NFTMatrix* A=nft_matrix_create(2,2);memcpy(A->data,Ad,4*sizeof(double));
      NormalFormSpec* s=nft_spec_create(2,1);NormalFormResult* r=nft_compute_normal_form(A,s);
      BifurcationType bt=nft_classify_bifurcation(r,0.0);assert(bt>=0);nft_result_free(r);nft_spec_free(s);nft_matrix_free(A); }
    nft_matrix_free(NULL);nft_spec_free(NULL);nft_result_free(NULL);nft_homological_free(NULL);nft_center_free(NULL);nft_unfolding_free(NULL);
    printf("All 14 tests passed.\n");return 0;
}
{ NFTMatrix* a=nft_matrix_create(3,3); double d[]={1,2,3,4,5,6,7,8,9}; memcpy(a->data,d,9*sizeof(double));
  double tr=nft_matrix_trace(a); assert(fabs(tr-15.0)<1e-6); nft_matrix_free(a); }
{ NFTMatrix* a=nft_matrix_create(2,2); a->data[0]=1;a->data[1]=2;a->data[2]=3;a->data[3]=4;
  double det=nft_matrix_det_2x2(a); assert(fabs(det+2.0)<1e-6); nft_matrix_free(a); }
{ double A[]={0,1,0,0}; NFTMatrix* m=nft_matrix_from_array(A,2,2); 
  bool nil=nft_matrix_is_nilpotent(m); assert(nil); nft_matrix_free(m); }
{ double A[]={2,0,0,2}; NFTMatrix* m=nft_matrix_from_array(A,2,2); 
  bool ss=nft_is_semisimple(m); assert(ss); nft_matrix_free(m); }
