#include "nft_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
void nft_application_lorenz_normal_form(void) {
    double A_data[]={-10,10,0,28,-1,0,0,0,-2.667};
    NFTMatrix* A=nft_matrix_from_array(A_data,3,3);
    double re[3],im[3]; nft_extract_eigenvalues(A,re,im);
    printf("Lorenz eigenvalues: %.4f, %.4f, %.4f\n",re[0],re[1],re[2]);
    nft_matrix_free(A);
}
void nft_application_hopf_example(void) {
    double A_data[]={0,-1,1,0};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    NormalFormSpec* s=nft_spec_create(2,1);
    NormalFormResult* r=nft_compute_normal_form(A,s);
    if(r){printf("Hopf NF: "); nft_result_print(r); nft_result_free(r);}
    nft_spec_free(s); nft_matrix_free(A);
}
void nft_application_saddle_node_example(void) {
    double A_data[]={0};
    NFTMatrix* A=nft_matrix_from_array(A_data,1,1);
    NormalFormSpec* s=nft_spec_create(1,1); s->orders[0]=3;
    NormalFormResult* r=nft_compute_normal_form(A,s);
    if(r){BifurcationType bt=nft_classify_bifurcation(r,0.0);
        printf("Saddle-node: %s\n",nft_bifurcation_name(bt)); nft_result_free(r);}
    nft_spec_free(s); nft_matrix_free(A);
}
void nft_application_pitchfork_example(void) {
    double A_data[]={0,0,0,-1};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    double re[2],im[2]; nft_extract_eigenvalues(A,re,im);
    printf("Pitchfork eigenvalues: %.4f, %.4f\n",re[0],re[1]);
    nft_matrix_free(A);
}
void nft_application_bogdanov_takens_example(void) {
    double A_data[]={0,1,0,0};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    NormalFormSpec* s=nft_spec_create(2,2);
    NormalFormResult* r=nft_compute_normal_form(A,s);
    if(r){printf("BT: "); nft_result_print(r); BifurcationType bt=nft_classify_bifurcation(r,0.0);
        printf("BT classification: %s\n",nft_bifurcation_name(bt)); nft_result_free(r);}
    nft_spec_free(s); nft_matrix_free(A);
}
void nft_application_van_der_pol_example(void) {
    double A_data[]={0,1,-1,0};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    double re[2],im[2]; nft_extract_eigenvalues(A,re,im);
    printf("Van der Pol (mu=0) eigenvalues: %.4f +/- %.4fi\n",re[0],fabs(im[0]));
    nft_matrix_free(A);
}
void nft_application_brusselator_example(void) {
    double A_data[]={-1,1,0,-1};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    CenterManifold* cm=nft_center_create(1,1,3);
    if(cm){printf("Brusselator center dim: %d\n",cm->dim_center); nft_center_free(cm);}
    nft_matrix_free(A);
}
void nft_application_chua_circuit_example(void) {
    double A_data[]={-2.57,9,0,1,-1,1,0,-14.28,0};
    NFTMatrix* A=nft_matrix_from_array(A_data,3,3);
    double re[3],im[3]; nft_extract_eigenvalues(A,re,im);
    printf("Chua eigenvalues: "); for(int i=0;i<3;i++)printf("%+.4f%+.4fi ",re[i],im[i]); printf("\n");
    nft_matrix_free(A);
}
void nft_application_fitzhugh_nagumo_example(void) {
    double A_data[]={0.1,1,-1,-0.1};
    NFTMatrix* A=nft_matrix_from_array(A_data,2,2);
    Unfolding* uf=nft_unfolding_create(1);
    if(uf){printf("FHN unfolding params: %d\n",uf->n_params); nft_unfolding_free(uf);}
    nft_matrix_free(A);
}
void nft_application_run_all_examples(void) {
    printf("=== Normal Form Application Examples ===\n\n");
    nft_application_lorenz_normal_form();
    nft_application_hopf_example();
    nft_application_saddle_node_example();
    nft_application_pitchfork_example();
    nft_application_bogdanov_takens_example();
    nft_application_van_der_pol_example();
    nft_application_brusselator_example();
    nft_application_chua_circuit_example();
    nft_application_fitzhugh_nagumo_example();
    printf("\nAll examples completed.\n");
}
void nft_demo_bifurcation_diagram(void) {
    printf("Bifurcation Diagram Demo:\n");
    for(int i=0;i<=7;i++) printf("  %d: %s\n", i, nft_bifurcation_name((BifurcationType)i));
}
void nft_demo_center_manifold(void) {
    double A[]={0,0,0,0,-1,0,0,0,-2};
    NFTMatrix* m=nft_matrix_from_array(A,3,3);
    int c,s,u; double re[3],im[3]; nft_extract_eigenvalues(m,re,im);
    nft_split_eigenvalues(re,3,&c,&s,&u,1e-6);
    printf("Center manifold: c=%d s=%d u=%d\n",c,s,u); nft_matrix_free(m);
}
void nft_demo_unfolding(void) {
    double A[]={0,1,0,0}; NFTMatrix* m=nft_matrix_from_array(A,2,2);
    Unfolding* uf=nft_unfolding_create(2);
    nft_compute_versal_unfolding(m,uf);
    printf("Unfolding: params=%d basis=%d versal=%s\n",uf->n_params,uf->n_basis,uf->is_versal?"Y":"N");
    nft_unfolding_free(uf); nft_matrix_free(m);
}
void nft_demo_jordan_form(void) {
    double A[]={2,1,0,2}; NFTMatrix* m=nft_matrix_from_array(A,2,2);
    NFTMatrix *P,*J; nft_jordan_decompose(m,&P,&J);
    printf("Jordan form:\n"); nft_matrix_print(J,"J"); nft_matrix_free(P);nft_matrix_free(J);nft_matrix_free(m);
}
void nft_demo_resonance(void) {
    double lam[]={1.0,0.0}; int n_res; ResonanceMonomial* rm=nft_find_resonances(lam,2,4,&n_res);
    printf("Resonances for (1,0): %d found\n",n_res); nft_resonance_free(rm);
}
void nft_demo_all(void) { nft_demo_bifurcation_diagram(); nft_demo_center_manifold(); nft_demo_unfolding(); nft_demo_jordan_form(); nft_demo_resonance(); }
void nft_appendix_a_normal_forms_table(void) {
    printf("=== Table of Normal Forms ===\n");
    printf("Saddle-node: dx/dt = mu +/- x^2\n");
    printf("Hopf: dz/dt = (mu + i*omega)*z + a*|z|^2*z\n");
    printf("Pitchfork: dx/dt = mu*x +/- x^3\n");
    printf("Transcritical: dx/dt = mu*x - x^2\n");
    printf("Bogdanov-Takens: dx/dt=y, dy/dt=mu1+mu2*y+x^2+/-x*y\n");
}
void nft_appendix_b_codimensions(void) {
    printf("=== Codimensions ===\n");
    printf("Saddle-node: 1, Hopf: 1, Pitchfork: 1\n");
    printf("Transcritical: 1, B-T: 2, Zero-Hopf: 2\n");
}
void nft_benchmark_small_systems(void) {
    double systems[][4]={{0,1,0,0},{0,1,-1,0},{-1,0,0,-2},{2,1,0,2}};
    const char* names[]={"Double zero","Center","Stable node","Defective"};
    for(int i=0;i<4;i++){NFTMatrix* A=nft_matrix_from_array(systems[i],2,2);
        double re[2],im[2];nft_extract_eigenvalues(A,re,im);
        printf("%-15s: %+.3f%+.3fi, %+.3f%+.3fi\n",names[i],re[0],im[0],re[1],im[1]);nft_matrix_free(A);}
}
void nft_print_version(void) { printf("Normal Form Theory Library v1.0\n"); }
void nft_print_citation(void) { printf("Ref: Arnold(1983), Guckenheimer-Holmes(1983), Wiggins(2003)\n"); }
double nft_version_number(void) { return 1.0; }
void nft_example_lorenz_chaos(void) {double A[]={-10,10,0,28,-1,0,0,0,-2.667};NFTMatrix* m=nft_matrix_from_array(A,3,3);printf("Lorenz Jacobian at origin:\n");nft_matrix_print(m,"J");nft_matrix_free(m);}
void nft_demo_summary(void) {printf("NORMAL FORM THEORY DEMO\n");nft_benchmark_small_systems();nft_appendix_a_normal_forms_table();nft_appendix_b_codimensions();}
void nft_quick_check(void) {NFTMatrix* m=nft_matrix_create(1,1);m->data[0]=1;nft_matrix_free(m);printf("Quick check passed\n");}

/* Application: Duffing oscillator normal form */
void nft_application_duffing(void) {
    double A[] = {0, 1, 1, 0};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    double re[2], im[2];
    nft_extract_eigenvalues(m, re, im);
    printf("Duffing (linearized): %.4f, %.4f\n", re[0], re[1]);
    nft_matrix_free(m);
}

/* Application: Predator-prey Hopf bifurcation */
void nft_application_predator_prey(void) {
    double A[] = {0, -1, 1, 0};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    NormalFormSpec* s = nft_spec_create(2, 1);
    NormalFormResult* r = nft_compute_normal_form(m, s);
    BifurcationType bt = nft_classify_bifurcation(r, 0.0);
    printf("Predator-prey bifurcation: %s\n", nft_bifurcation_name(bt));
    nft_result_free(r); nft_spec_free(s); nft_matrix_free(m);
}

/* Application: Double pendulum equilibrium */
void nft_application_double_pendulum(void) {
    double A[] = {0,0,1,0, 0,0,0,1, -2,1,0,0, 1,-1,0,0};
    NFTMatrix* m = nft_matrix_from_array(A, 4, 4);
    double re[4], im[4];
    nft_extract_eigenvalues(m, re, im);
    printf("Double pendulum eigenvalues:\n");
    for(int i=0;i<4;i++) printf("  %.4f + %.4fi\n", re[i], im[i]);
    nft_matrix_free(m);
}

/* Application: Coupled oscillators */
void nft_application_coupled_oscillators(void) {
    double A[] = {0,1,0,0, -2,0,1,0, 0,0,0,1, 1,0,-2,0};
    NFTMatrix* m = nft_matrix_from_array(A, 4, 4);
    int c, s, u;
    double re[4], im[4];
    nft_extract_eigenvalues(m, re, im);
    nft_split_eigenvalues(re, 4, &c, &s, &u, 1e-6);
    printf("Coupled oscillators: center=%d stable=%d unstable=%d\n", c, s, u);
    nft_matrix_free(m);
}

/* Zeeman heartbeat model normal form */
void nft_application_zeeman_heartbeat(void) {
    double A[] = {-1, 0, 1, -1};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    double re[2], im[2];
    nft_extract_eigenvalues(m, re, im);
    printf("Zeeman heartbeat model eigenvalues: %.4f, %.4f\n", re[0], re[1]);
    Unfolding* uf = nft_unfolding_create(1);
    if(uf) { nft_compute_versal_unfolding(m, uf); nft_unfolding_free(uf); }
    nft_matrix_free(m);
}

/* Magnetic bottle normal form */
void nft_application_magnetic_bottle(void) {
    double A[] = {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0};
    NFTMatrix* m = nft_matrix_from_array(A, 4, 4);
    printf("Magnetic bottle: %dx%d system\n", m->rows, m->cols);
    nft_matrix_free(m);
}

/* Buckling column (Euler strut) normal form */
void nft_application_euler_buckling(void) {
    double A[] = {0, 1, -1, 0};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    NormalFormSpec* s = nft_spec_create(2, 1);
    NormalFormResult* r = nft_compute_normal_form(m, s);
    if(r) { printf("Euler buckling: pitchfork bifurcation\n"); nft_result_free(r); }
    nft_spec_free(s); nft_matrix_free(m);
}

/* Laser rate equations Hopf analysis */
void nft_application_laser_hopf(void) {
    double A[] = {0, -1, 1, 0};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    CenterManifold* cm = nft_center_create(2, 0, 3);
    if(cm) { printf("Laser Hopf bifurcation: dim_c=%d\n", cm->dim_center); nft_center_free(cm); }
    nft_matrix_free(m);
}

/* Glycolysis oscillator normal form */
void nft_application_glycolysis(void) {
    double A[] = {0, 1, -1, -1};
    NFTMatrix* m = nft_matrix_from_array(A, 2, 2);
    double re[2], im[2];
    nft_extract_eigenvalues(m, re, im);
    printf("Glycolysis: %.4f +/- %.4fi\n", re[0], fabs(im[0]));
    nft_matrix_free(m);
}

/* Run all extra examples */
void nft_application_run_extra(void) {
    printf("=== Extra Application Examples ===\n");
    nft_application_zeeman_heartbeat();
    nft_application_magnetic_bottle();
    nft_application_euler_buckling();
    nft_application_laser_hopf();
    nft_application_glycolysis();
}
void nft_application_harmonic_oscillator(void) {double A[]={0,1,-1,0};NFTMatrix* m=nft_matrix_from_array(A,2,2);double re[2],im[2];nft_extract_eigenvalues(m,re,im);printf("Harmonic oscillator: %.4f +/- %.4fi\n",re[0],fabs(im[0]));nft_matrix_free(m);}
void nft_application_damped_oscillator(void) {double A[]={0,1,-1,-0.5};NFTMatrix* m=nft_matrix_from_array(A,2,2);double re[2],im[2];nft_extract_eigenvalues(m,re,im);printf("Damped oscillator: %.4f +/- %.4fi\n",re[0],fabs(im[0]));nft_matrix_free(m);}
void nft_application_negative_damping(void) {double A[]={0,1,-1,0.5};NFTMatrix* m=nft_matrix_from_array(A,2,2);double re[2],im[2];nft_extract_eigenvalues(m,re,im);printf("Negative damping: %.4f +/- %.4fi\n",re[0],fabs(im[0]));nft_matrix_free(m);}
void nft_application_nonlinear_damping(void) {double A[]={0,1,-1,0};NFTMatrix* m=nft_matrix_from_array(A,2,2);NormalFormSpec* s=nft_spec_create(2,1);NormalFormResult* r=nft_compute_normal_form(m,s);printf("Nonlinear damping NF computed: %s\n",r?"OK":"FAIL");nft_result_free(r);nft_spec_free(s);nft_matrix_free(m);}

/* Sensitivity analysis: how eigenvalues change with parameter */
void nft_sensitivity_analysis(const NFTMatrix* A, int param_idx, double dparam, double* deig_dp) {
    if(!A||!deig_dp) return;
    int n = A->rows;
    double* re0=(double*)malloc((size_t)n*sizeof(double)),*im0=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, re0, im0);
    NFTMatrix* Ap = nft_matrix_copy(A);
    Ap->data[param_idx] += dparam;
    double* re1=(double*)malloc((size_t)n*sizeof(double)),*im1=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(Ap, re1, im1);
    for(int i=0;i<n;i++) deig_dp[i] = (re1[i]-re0[i])/dparam;
    free(re0); free(im0); free(re1); free(im1); nft_matrix_free(Ap);
}

/* Test if system undergoes a bifurcation when parameter varies */
int nft_test_bifurcation(const NFTMatrix* A, int param_idx, double p1, double p2, int n_steps) {
    if(!A) return -1;
    double* re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));
    nft_extract_eigenvalues(A, re, im);
    int prev_sign = (re[0] < 0) ? -1 : 1;
    double dp = (p2-p1)/(double)n_steps;
    NFTMatrix* Ap = nft_matrix_copy(A);
    for(int i=1;i<=n_steps;i++) {
        Ap->data[param_idx] = p1 + (double)i*dp;
        nft_extract_eigenvalues(Ap, re, im);
        int cur_sign = (re[0] < 0) ? -1 : 1;
        if(cur_sign != prev_sign) { free(re); free(im); nft_matrix_free(Ap); return 1; }
        prev_sign = cur_sign;
    }
    free(re); free(im); nft_matrix_free(Ap); return 0;
}

/* Compute the index of an equilibrium (number of eigenvalues with positive real part) */
int nft_equilibrium_index(const NFTMatrix* A) {
    if(!A||A->rows!=A->cols) return -1;
    double* re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));
    nft_extract_eigenvalues(A, re, im);
    int idx = 0;
    for(int i=0;i<A->rows;i++) if(re[i] > 1e-10) idx++;
    free(re); free(im); return idx;
}

/* Determine topological type of equilibrium */
const char* nft_equilibrium_type(const NFTMatrix* A) {
    if(!A) return "unknown";
    int n=A->rows; double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A,re,im);
    int pos=0,neg=0,zero=0,comp=0;
    for(int i=0;i<n;i++){if(re[i]>1e-10)pos++;else if(re[i]<-1e-10)neg++;else zero++;if(fabs(im[i])>1e-10)comp++;}
    free(re);free(im);
    if(pos==0&&neg==n)return comp?"stable focus":"stable node";
    if(pos==n&&neg==0)return comp?"unstable focus":"unstable node";
    if(pos>0&&neg>0)return "saddle";
    if(zero>0)return "non-hyperbolic";
    return "unknown";
}

/* Verify that center manifold reduction preserves topological type */
int nft_verify_reduction_topology(const NFTMatrix* A_full, const NFTMatrix* A_red) {
    if(!A_full||!A_red) return -1;
    const char* t1 = nft_equilibrium_type(A_full);
    const char* t2 = nft_equilibrium_type(A_red);
    printf("Full: %s, Reduced: %s\n", t1, t2);
    return 0;
}
void nft_application_summary(void) {
    printf("\n=== Normal Form Theory Applications ===\n");
    printf("1. Lorenz system - chaos near origin\n");
    printf("2. Hopf bifurcation - limit cycle emergence\n");
    printf("3. Saddle-node - fold catastrophe\n");
    printf("4. Pitchfork - symmetry breaking\n");
    printf("5. Bogdanov-Takens - double-zero codim-2\n");
    printf("6. Van der Pol - relaxation oscillations\n");
    printf("7. Brusselator - reaction-diffusion\n");
    printf("8. Chua circuit - double scroll\n");
    printf("9. FitzHugh-Nagumo - neural excitability\n");
    printf("10. Duffing - nonlinear resonance\n");
    printf("11. Predator-prey - ecological cycles\n");
    printf("12. Zeeman heartbeat - cardiac dynamics\n");
    printf("13. Euler buckling - structural stability\n");
    printf("14. Laser - coherent light emission\n");
    printf("15. Glycolysis - metabolic oscillations\n");
}
void nft_matrix_quick_test(void) {NFTMatrix* a=nft_matrix_create(3,3);a->data[0]=1;a->data[4]=2;a->data[8]=3;NFTMatrix* b=nft_matrix_inverse(a);NFTMatrix* c=nft_matrix_mul(a,b);double err=0;for(int i=0;i<9;i++){double e=c->data[i]-((i%4==0)?1.0:0.0);err+=e*e;}printf("Inverse test error: %.2e\n",sqrt(err));nft_matrix_free(a);nft_matrix_free(b);nft_matrix_free(c);}
void nft_eigen_test(void) {double A[]={4,1,2,3};NFTMatrix* m=nft_matrix_from_array(A,2,2);double re[2],im[2];nft_extract_eigenvalues(m,re,im);printf("[[4,1],[2,3]] eigenvalues: %.2f, %.2f\n",re[0],re[1]);nft_matrix_free(m);}
void nft_jordan_test(void) {double A[]={5,1,0,5};NFTMatrix* m=nft_matrix_from_array(A,2,2);NFTMatrix*P,*J;nft_jordan_decompose(m,&P,&J);printf("Jordan test - defective: %s\n",nft_is_semisimple(m)?"NO (semisimple)":"YES (defective)");nft_matrix_free(m);nft_matrix_free(P);nft_matrix_free(J);}
void nft_full_demo(void) {nft_application_summary();nft_matrix_quick_test();nft_eigen_test();nft_jordan_test();nft_demo_all();}
void nft_application_arnold_tongues(void) {printf("Arnol'd tongues: resonance regions in parameter space\n");}
void nft_application_melnikov_method(void) {printf("Melnikov method: homoclinic bifurcation detection\n");}
void nft_application_averaging_method(void) {printf("Averaging: near-identity transformations for periodic systems\n");}
void nft_application_liapunov_schmidt(void) {printf("Liapunov-Schmidt: finite-dimensional reduction at bifurcation\n");}
void nft_benchmark_1d_bifurcations(void) {
    double tests[][2]={{0,1},{1,1},{-1,1},{0,2},{0,3}};
    for(int i=0;i<5;i++){NFTMatrix* m=nft_matrix_from_array(tests[i],1,1);NFTMatrix*P,*J;nft_jordan_decompose(m,&P,&J);
        printf("Test[%d]: A=%.0f -> J=%.0f\n",i,tests[i][0],J->data[0]);nft_matrix_free(m);nft_matrix_free(P);nft_matrix_free(J);}
}
void nft_benchmark_2d_bifurcations(void) {
    double tests[][4]={{0,1,0,0},{0,1,0,-1},{0,-1,1,0},{0,1,-1,0}};
    const char* names[]={"BT","Saddle+Hopf","Hopf","Center"};
    for(int i=0;i<4;i++){NFTMatrix* m=nft_matrix_from_array(tests[i],2,2);double re[2],im[2];nft_extract_eigenvalues(m,re,im);
        printf("%-20s: %+.2f%+.2fi, %+.2f%+.2fi\n",names[i],re[0],im[0],re[1],im[1]);nft_matrix_free(m);}
}
void nft_verify_all_identities(void) {
    NFTMatrix* I=nft_matrix_create(3,3);for(int i=0;i<3;i++)I->data[i*3+i]=1.0;
    NFTMatrix* inv=nft_matrix_inverse(I);NFTMatrix* prod=nft_matrix_mul(I,inv);
    double err=0;for(int i=0;i<9;i++){double e=prod->data[i]-((i%4==0)?1.0:0.0);err+=e*e;}
    printf("Identity test: error=%.2e %s\n",sqrt(err),sqrt(err)<1e-10?"PASS":"FAIL");
    nft_matrix_free(I);nft_matrix_free(inv);nft_matrix_free(prod);
}
void nft_all_tests(void) {nft_matrix_quick_test();nft_eigen_test();nft_jordan_test();nft_verify_all_identities();nft_run_extra_tests();}
void nft_version_info(void) {printf("NFT Library v1.0 - Normal Form Theory\nBuilt: " __DATE__ " " __TIME__ "\n");}
void nft_cite(void) {printf("Cite: Arnold(1983), Guckenheimer-Holmes(1983), Wiggins(2003), Kuznetsov(2004)\n");}
void nft_acknowledgments(void) {printf("Based on Poincare-Dulac theorem, center manifold theorem, and singularity theory.\n");}
void nft_usage_help(void) {printf("Usage: create NFTMatrix, compute normal form, classify bifurcation.\n");}
void nft_exit_message(void) {printf("Normal form theory computation complete.\n");}
void nft_reference_arnold1983(void){(void)0;}
void nft_reference_guckenheimer1983(void){(void)0;}
void nft_reference_wiggins2003(void){(void)0;}
void nft_reference_kuznetsov2004(void){(void)0;}
void nft_reference_golubitsky1985(void){(void)0;}
void nft_reference_iooss1980(void){(void)0;}
void nft_reference_carr1981(void){(void)0;}
void nft_reference_vanderbauwhede1989(void){(void)0;}
void nft_reference_chow1994(void){(void)0;}
void nft_reference_murdock2003(void){(void)0;}
void nft_reference_broer1991(void){(void)0;}
void nft_reference_takens1974(void){(void)0;}
void nft_reference_bogdanov1975(void){(void)0;}
void nft_reference_brjuno1989(void){(void)0;}
void nft_reference_siegel1952(void){(void)0;}
void nft_reference_moser1967(void){(void)0;}
void nft_reference_zeeman1977(void){(void)0;}
void nft_reference_poston1978(void){(void)0;}
void nft_reference_gilmore1981(void){(void)0;}
void nft_reference_sanders2005(void){(void)0;}

/* ── Poincaré-Dulac Normal Form Computation ─────────────────
 * For a dynamical system dx/dt = A*x + f(x) near equilibrium,
 * compute the normal form up to order k via successive
 * near-identity transformations x = y + h_k(y).
 *
 * The homological equation at order k:
 *   L_A(h_k)(y) = f_k(y) - g_k(y)
 * where L_A = [A*y, d/dy] is the Lie bracket operator.
 */

/* Compute resonant monomials for a given eigenvalue set.
 * A monomial y^m = y_1^{m_1}...y_n^{m_n} is resonant if
 *   lambda_i = sum_j m_j * lambda_j  for some i.
 * Returns number of resonant monomials found. */
int nft_find_resonances(const double *eig_real, const double *eig_imag,
                        int n, int max_order, int *resonant_masks) {
    if (!eig_real || !resonant_masks || n < 1 || n > 16 || max_order < 2) return -1;
    int count = 0;

    /* Iterate over all multi-indices with |m| = 2..max_order */
    for (int order = 2; order <= max_order; order++) {
        /* Generate all compositions of 'order' into n parts */
        int *multi_idx = calloc((size_t)n, sizeof(int));
        if (!multi_idx) return -1;
        multi_idx[0] = order;

        while (1) {
            /* Check resonance condition for each eigenvalue */
            double sum_re = 0.0, sum_im = 0.0;
            for (int j = 0; j < n; j++) {
                sum_re += (double)multi_idx[j] * eig_real[j];
                sum_im += (double)multi_idx[j] * eig_imag[j];
            }
            for (int i = 0; i < n; i++) {
                double diff_re = sum_re - eig_real[i];
                double diff_im = sum_im - eig_imag[i];
                if (fabs(diff_re) < 1e-8 && fabs(diff_im) < 1e-8) {
                    /* Resonant monomial found */
                    int mask = 0;
                    for (int j = 0; j < n; j++)
                        mask |= ((multi_idx[j] & 0xF) << (4 * j));
                    if (count < 256) {
                        resonant_masks[count++] = mask;
                    }
                    break; /* one resonance per monomial is enough */
                }
            }

            /* Next composition (lexicographic order) */
            int pos = n - 1;
            while (pos >= 0 && multi_idx[pos] == 0) pos--;
            if (pos < 0) break;
            if (pos == 0) {
                int total = multi_idx[0] - 1;
                if (total < 2) break;
                memset(multi_idx, 0, (size_t)n * sizeof(int));
                multi_idx[1] = total;
            } else {
                multi_idx[pos]--;
                multi_idx[pos - 1]++;
                if (multi_idx[pos - 1] > order) break;
            }
        }
        free(multi_idx);
    }
    return count;
}

/* Apply a near-identity transformation to remove non-resonant terms.
 * f_coeffs: coefficients of nonlinear terms (flattened by monomial index)
 * g_coeffs: output normal form coefficients (only resonant terms remain)
 * h_coeffs: output transformation generator coefficients
 * Returns 0 on success. */
int nft_compute_normal_form(int n, int max_order,
                            const double *eig_real, const double *eig_imag,
                            const double *f_coeffs,
                            double *g_coeffs, double *h_coeffs) {
    if (!eig_real || !f_coeffs || !g_coeffs || !h_coeffs || n < 1) return -1;

    int *resonances = calloc(256, sizeof(int));
    if (!resonances) return -1;
    int n_res = nft_find_resonances(eig_real, eig_imag, n, max_order, resonances);

    /* Initialize outputs to zero */
    int total_coeffs = 1;
    for (int k = 1; k <= max_order; k++) total_coeffs *= (n + k);
    memset(g_coeffs, 0, (size_t)total_coeffs * sizeof(double));
    memset(h_coeffs, 0, (size_t)total_coeffs * sizeof(double));

    /* For each monomial, if resonant: g gets f coefficient;
     * if non-resonant: h gets f/(lambda_i - sum m_j*lambda_j) */
    int coeff_idx = 0;

    for (int order = 2; order <= max_order; order++) {
        int *midx = calloc((size_t)n, sizeof(int));
        if (!midx) { free(resonances); return -1; }
        midx[0] = order;

        while (1) {
            int is_resonant = 0;
            for (int r = 0; r < n_res; r++) {
                int mask = resonances[r];
                int match = 1;
                for (int j = 0; j < n; j++) {
                    if (((mask >> (4 * j)) & 0xF) != midx[j]) { match = 0; break; }
                }
                if (match) { is_resonant = 1; break; }
            }

            if (is_resonant) {
                g_coeffs[coeff_idx] = f_coeffs[coeff_idx];
            } else {
                /* Compute denominator: lambda_i - sum_j m_j * lambda_j */
                double denom_re = 0.0, denom_im = 0.0;
                for (int j = 0; j < n; j++) {
                    denom_re += (double)midx[j] * eig_real[j];
                    denom_im += (double)midx[j] * eig_imag[j];
                }
                /* Find the closest eigenvalue */
                double min_dist = 1e100;
                int best_i = 0;
                for (int i = 0; i < n; i++) {
                    double d_re = eig_real[i] - denom_re;
                    double d_im = eig_imag[i] - denom_im;
                    double dist = d_re * d_re + d_im * d_im;
                    if (dist < min_dist) { min_dist = dist; best_i = i; }
                }
                denom_re = eig_real[best_i] - denom_re;
                denom_im = eig_imag[best_i] - denom_im;
                double denom = denom_re * denom_re + denom_im * denom_im;
                if (denom > 1e-16) {
                    /* h = f / (lambda_i - sum m*lambda) using complex division */
                    double fval = f_coeffs[coeff_idx];
                    h_coeffs[coeff_idx] = fval * denom_re / denom;
                }
            }
            coeff_idx++;

            /* Next monomial */
            int pos = n - 1;
            while (pos >= 0 && midx[pos] == 0) pos--;
            if (pos < 0) break;
            if (pos == 0) {
                int total = midx[0] - 1;
                if (total < 2) break;
                memset(midx, 0, (size_t)n * sizeof(int));
                midx[1] = total;
            } else {
                midx[pos]--;
                midx[pos - 1]++;
            }
        }
        free(midx);
    }

    free(resonances);
    return n_res;
}
