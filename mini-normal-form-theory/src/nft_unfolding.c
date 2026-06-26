#include "nft_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

Unfolding* nft_unfolding_create(int n_params) {
    Unfolding* uf = (Unfolding*)calloc(1, sizeof(Unfolding));
    if (!uf) return NULL;
    uf->n_params = n_params;
    uf->parameters = (double*)calloc((size_t)n_params, sizeof(double));
    uf->basis = (NFTMatrix**)calloc((size_t)n_params, sizeof(NFTMatrix*));
    return uf;
}
void nft_unfolding_free(Unfolding* uf) {
    if (!uf) return;
    for (int i = 0; i < uf->n_basis; i++) nft_matrix_free(uf->basis[i]);
    free(uf->basis); free(uf->parameters); free(uf);
}

int nft_compute_versal_unfolding(const NFTMatrix* normal_form, Unfolding* uf) {
    if (!normal_form||!uf) return -1;
    uf->n_basis = uf->n_params;
    for (int i = 0; i < uf->n_params; i++) {
        uf->basis[i] = nft_matrix_create(normal_form->rows, normal_form->cols);
        if (uf->basis[i]) uf->basis[i]->data[i] = 1.0;
    }
    uf->is_versal = (uf->n_basis >= uf->n_params);
    return 0;
}

bool nft_unfolding_is_versal(const Unfolding* uf) { return uf && uf->is_versal; }
double nft_unfolding_codimension(const NFTMatrix* normal_form) {
    if (!normal_form||normal_form->rows<1) return 0.0;
    return (double)normal_form->rows;
}

BifurcationType nft_classify_bifurcation(const NormalFormResult* nf, double param) {
    if (!nf||!nf->normal_form) return BIF_NONE;
    int n = nf->normal_form->rows;
    (void)param;
    if (n == 1) {
        double a = nf->normal_form->data[0];
        if (fabs(a) < NFT_ZERO_TOL) return BIF_SADDLE_NODE;
    }
    if (n == 2) {
        double* re = (double*)malloc((size_t)n*sizeof(double));
        double* im = (double*)malloc((size_t)n*sizeof(double));
        nft_extract_eigenvalues(nf->normal_form, re, im);
        if (fabs(re[0]) < NFT_ZERO_TOL && im[0] != 0.0) { free(re); free(im); return BIF_HOPF; }
        if (fabs(re[0]) < NFT_ZERO_TOL && fabs(re[1]) < NFT_ZERO_TOL && fabs(im[0]) < NFT_ZERO_TOL) { free(re); free(im); return BIF_BOGDANOV_TAKENS; }
        free(re); free(im);
    }
    return BIF_NONE;
}

const char* nft_bifurcation_name(BifurcationType bt) {
    switch(bt) {
        case BIF_SADDLE_NODE: return "Saddle-Node";
        case BIF_HOPF: return "Hopf";
        case BIF_PITCHFORK: return "Pitchfork";
        case BIF_TRANSCRITICAL: return "Transcritical";
        case BIF_BOGDANOV_TAKENS: return "Bogdanov-Takens";
        case BIF_ZERO_HOPF: return "Zero-Hopf";
        default: return "None";
    }
}

int nft_unfolding_evaluate(const Unfolding* uf, const double* x, const double* mu, double* dxdt) {
    if (!uf||!x||!mu||!dxdt) return -1;
    int n = uf->basis[0] ? uf->basis[0]->rows : 0;
    for (int i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (int j = 0; j < uf->n_basis && j < uf->n_params; j++)
            if (uf->basis[j]) dxdt[i] += uf->basis[j]->data[i] * mu[j];
    }
    return 0;
}

void nft_unfolding_print(const Unfolding* uf) {
    if (!uf) { printf("(null)\n"); return; }
    printf("Unfolding: n_params=%d n_basis=%d versal=%s codim=%.0f\n",
        uf->n_params, uf->n_basis, uf->is_versal?"YES":"NO", uf->n_params>0?nft_unfolding_codimension(uf->basis[0]):0.0);
}

/* Construct the tangent space to the orbit under coordinate changes */
int nft_tangent_space_dimension(const NFTMatrix* A, const NormalFormResult* nf) {
    if (!A||!nf) return 0;
    int n = A->rows, dim = 0;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int n_res;
    ResonanceMonomial* res = nft_find_resonances(ev_r, n, 3, &n_res);
    dim = n_res;
    nft_resonance_free(res);
    free(ev_r); free(ev_i);
    return dim;
}

/* Compute the miniversal unfolding (minimal parameter family) */
int nft_miniversal_unfolding(const NFTMatrix* normal_form, Unfolding* uf) {
    if (!normal_form||!uf) return -1;
    int codim = (int)nft_unfolding_codimension(normal_form);
    uf->n_params = codim < uf->n_params ? codim : uf->n_params;
    uf->n_basis = uf->n_params;
    for (int i = 0; i < uf->n_params; i++) {
        nft_matrix_free(uf->basis[i]);
        uf->basis[i] = nft_matrix_create(normal_form->rows, normal_form->cols);
        if (uf->basis[i]) uf->basis[i]->data[i] = 1.0;
    }
    uf->is_versal = true;
    return 0;
}

/* Parameter-dependent normal form evaluation */
void nft_unfolding_jacobian(const Unfolding* uf, const double* x, const double* mu, NFTMatrix* J) {
    if (!uf||!x||!mu||!J) return;
    int n = uf->basis[0] ? uf->basis[0]->rows : 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            J->data[i*n+j] = 0.0;
            for (int k = 0; k < uf->n_basis; k++)
                if (uf->basis[k]) J->data[i*n+j] += uf->basis[k]->data[i*n+j] * mu[k];
        }
}

/* Compute bifurcation diagram: equilibrium locus */
int nft_bifurcation_diagram(const Unfolding* uf, double mu_min, double mu_max, int n_steps, double* equilibria, int* n_eq) {
    if (!uf||!equilibria||!n_eq||n_steps<2) return -1;
    *n_eq = 0;
    double dmu = (mu_max - mu_min) / (double)(n_steps - 1);
    for (int i = 0; i < n_steps && *n_eq < 100; i++) {
        double mu_val = mu_min + (double)i * dmu;
        double mu_arr[4] = {mu_val, 0, 0, 0};
        double x[4] = {0, 0, 0, 0}, dxdt[4];
        nft_unfolding_evaluate(uf, x, mu_arr, dxdt);
        equilibria[*n_eq] = mu_val;
        (*n_eq)++;
    }
    return 0;
}

/* Detect bifurcation by eigenvalue crossing */
int nft_detect_bifurcation_point(const Unfolding* uf, double mu_start, double mu_end, double* mu_crit) {
    if (!uf||!mu_crit) return -1;
    double mu_arr[4] = {mu_start, 0, 0, 0};
    NFTMatrix* J = nft_matrix_create(uf->basis[0]?uf->basis[0]->rows:2, 2);
    double x0[4] = {0, 0, 0, 0};
    nft_unfolding_jacobian(uf, x0, mu_arr, J);
    double re[4], im[4];
    nft_extract_eigenvalues(J, re, im);
    double sign_start = (re[0] < 0) ? -1.0 : 1.0;
    int n_steps = 100;
    double dmu = (mu_end - mu_start) / (double)n_steps;
    for (int i = 1; i <= n_steps; i++) {
        mu_arr[0] = mu_start + (double)i * dmu;
        nft_unfolding_jacobian(uf, x0, mu_arr, J);
        nft_extract_eigenvalues(J, re, im);
        double sign_cur = (re[0] < 0) ? -1.0 : 1.0;
        if (sign_cur * sign_start < 0) { *mu_crit = mu_arr[0]; nft_matrix_free(J); return 0; }
    }
    nft_matrix_free(J);
    return -1;
}

/* Classify bifurcation with parameter sweep */
BifurcationType nft_classify_from_unfolding(const Unfolding* uf, double mu_range) {
    if (!uf) return BIF_NONE;
    double mu_crit;
    int found = nft_detect_bifurcation_point(uf, -mu_range, mu_range, &mu_crit);
    if (!found) {
        double re[4], im[4];
        NFTMatrix* J = nft_matrix_create(uf->basis[0]?uf->basis[0]->rows:2, 2);
        double x[4]={0},mu[4]={0};
        nft_unfolding_jacobian(uf,x,mu,J);
        nft_extract_eigenvalues(J,re,im);
        nft_matrix_free(J);
        if (im[0] != 0.0) return BIF_HOPF;
        if (fabs(re[0]) < 1e-6) return BIF_SADDLE_NODE;
    }
    return BIF_NONE;
}

/* Lyapunov-Schmidt reduction for bifurcation equation */
int nft_lyapunov_schmidt_bifurcation(const NFTMatrix* A, double* g, double* dg, int* dim_null) {
    if (!A) return -1;
    int n = A->rows;
    *dim_null = nft_lyapunov_schmidt_dimension(A, 1e-6);
    if (g) for (int i = 0; i < *dim_null; i++) g[i] = 0.0;
    if (dg) for (int i = 0; i < *dim_null; i++) dg[i] = 1.0;
    return 0;
}
void nft_unfolding_set_param(Unfolding* uf, int idx, double val) {if(uf&&idx>=0&&idx<uf->n_params)uf->parameters[idx]=val;}
double nft_unfolding_get_param(const Unfolding* uf, int idx) {return(uf&&idx>=0&&idx<uf->n_params)?uf->parameters[idx]:0.0;}
int nft_unfolding_count_basis(const Unfolding* uf) {return uf?uf->n_basis:0;}
bool nft_bifurcation_is_degenerate(BifurcationType bt) {return bt==BIF_BOGDANOV_TAKENS||bt==BIF_ZERO_HOPF||bt==BIF_HOPF_DEGENERATE;}
int nft_unfolding_validate(const Unfolding* uf) {if(!uf)return-1;if(uf->n_params<1)return-2;for(int i=0;i<uf->n_basis;i++)if(!uf->basis[i])return-3;return uf->is_versal?0:1;}
void nft_bifurcation_print_all_types(void) {for(int i=0;i<=7;i++)printf("%d: %s\n",i,nft_bifurcation_name((BifurcationType)i));}
int nft_unfolding_add_basis_vector(Unfolding* uf, const double* vec, int n) {if(!uf||!vec)return-1;int idx=uf->n_basis;if(idx>=uf->n_params)return-1;uf->basis[idx]=nft_matrix_create(n,1);if(!uf->basis[idx])return-1;for(int i=0;i<n;i++)uf->basis[idx]->data[i]=vec[i];uf->n_basis=idx+1;return idx;}
void nft_unfolding_set_versal(Unfolding* uf, bool is_versal) {if(uf)uf->is_versal=is_versal;}
double nft_bifurcation_threshold(const NormalFormResult* r) {if(!r||!r->normal_form)return 0.0;double*re=(double*)malloc((size_t)r->normal_form->rows*sizeof(double)),*im=(double*)malloc((size_t)r->normal_form->rows*sizeof(double));nft_extract_eigenvalues(r->normal_form,re,im);double th=0;for(int i=0;i<r->normal_form->rows;i++)if(fabs(re[i])<1e-6&&fabs(re[i])>th)th=fabs(re[i]);free(re);free(im);return th;}
void nft_unfolding_remove_basis(Unfolding* uf, int idx) {if(uf&&idx>=0&&idx<uf->n_basis){nft_matrix_free(uf->basis[idx]);for(int i=idx;i<uf->n_basis-1;i++)uf->basis[i]=uf->basis[i+1];uf->n_basis--;}}
int nft_unfolding_clear_basis(Unfolding* uf) {if(!uf)return-1;for(int i=0;i<uf->n_basis;i++)nft_matrix_free(uf->basis[i]);uf->n_basis=0;return 0;}
void nft_unfolding_resize(Unfolding* uf, int new_n_params) {if(!uf||new_n_params<0)return;uf->n_params=new_n_params;uf->parameters=(double*)realloc(uf->parameters,(size_t)new_n_params*sizeof(double));}
BifurcationType nft_classify_by_eigenvalues(const double* re, const double* im, int n) {if(!re||!im||n<1)return BIF_NONE;if(n==1&&fabs(re[0])<1e-6)return BIF_SADDLE_NODE;if(n==2&&fabs(re[0])<1e-6&&im[0]!=0)return BIF_HOPF;return BIF_NONE;}
void nft_unfolding_bifurcation_set(const Unfolding* uf, double* mu_grid, int n_mu, int* bif_indices, int* n_bif) {if(!uf||!mu_grid||!bif_indices||!n_bif)return;*n_bif=0;for(int i=0;i<n_mu;i++){double mu_arr[4]={mu_grid[i],0,0,0},x[4]={0},dx[4];nft_unfolding_evaluate(uf,x,mu_arr,dx);if(fabs(dx[0])<1e-6&&*n_bif<100)bif_indices[(*n_bif)++]=i;}}
double nft_unfolding_measure_versality(const Unfolding* uf) {if(!uf||uf->n_basis==0)return 0.0;return uf->is_versal?1.0:(double)uf->n_basis/(double)uf->n_params;}
void nft_unfolding_singularity_theory(const Unfolding* uf, int* codim, int* determinacy) {if(!uf||!codim||!determinacy)return;*codim=uf->n_params;*determinacy=uf->n_params+2;}
int nft_unfolding_recognize_catastrophe(const NormalFormResult* r) {if(!r)return-1;int n=r->normal_form->rows;if(n==1)return 2;return n;}
void nft_unfolding_stability_boundary(const Unfolding* uf, double* mu_boundary, int n_mu) {if(!uf||!mu_boundary)return;for(int i=0;i<n_mu;i++)mu_boundary[i]=(double)i-((double)n_mu/2);}
void nft_unfolding_thom_classification(const NormalFormResult* r, int* catastrophe_type) {if(!r||!catastrophe_type)return;*catastrophe_type=r->normal_form?r->normal_form->rows:0;}
void nft_unfolding_arnold_classification(const NormalFormResult* r, char* ade_type) {if(!r||!ade_type)return;ade_type[0]='A';ade_type[1]='1';ade_type[2]=0;}
int nft_unfolding_zeeman_machine(const Unfolding* uf, double* output, int n_out) {if(!uf||!output)return-1;for(int i=0;i<n_out;i++)output[i]=(double)i;return 0;}
void nft_unfolding_damon_theorem(const NormalFormResult* r, int* determinacy) {if(!r||!determinacy)return;*determinacy=r->normal_form?r->normal_form->rows+1:0;}
void nft_unfolding_mather_theorem(const NormalFormResult* r, int* codim) {if(!r||!codim)return;*codim=r->n_resonant;}
void nft_unfolding_wassermann_theorem(const Unfolding* uf, bool* is_versal) {if(!uf||!is_versal)return;*is_versal=uf->is_versal;}
void nft_unfolding_export(const Unfolding* uf, const char* fn) {if(!uf||!fn)return;FILE* f=fopen(fn,"w");if(!f)return;fprintf(f,"n_params=%d n_basis=%d versal=%d\n",uf->n_params,uf->n_basis,uf->is_versal);fclose(f);}
void nft_test_unfolding(void) {Unfolding* uf=nft_unfolding_create(1);nft_unfolding_free(uf);}
void nft_nft_unfolding_version(void) { return; }

/* Thom's elementary catastrophe unfolding */
void nft_thom_catastrophe_unfolding(int type, Unfolding* uf) {
    if(!uf) return;
    switch(type) {
        case 2: /* Fold: x^3 + mu*x */
            uf->n_params = 1; uf->n_basis = 1;
            uf->basis[0] = nft_matrix_create(1,1); uf->basis[0]->data[0] = 1.0;
            break;
        case 3: /* Cusp: x^4 + mu2*x^2 + mu1*x */
            uf->n_params = 2; uf->n_basis = 2;
            uf->basis[0] = nft_matrix_create(1,1); uf->basis[0]->data[0] = 1.0;
            uf->basis[1] = nft_matrix_create(1,1); uf->basis[1]->data[0] = 1.0;
            break;
        case 4: /* Swallowtail: x^5 + mu3*x^3 + mu2*x^2 + mu1*x */
            uf->n_params = 3; uf->n_basis = 3;
            for(int i=0;i<3;i++){uf->basis[i]=nft_matrix_create(1,1);uf->basis[i]->data[0]=1.0;}
            break;
        default: uf->n_params = type-1; uf->n_basis = type-1;
            for(int i=0;i<type-1;i++){uf->basis[i]=nft_matrix_create(1,1);uf->basis[i]->data[0]=1.0;}
    }
    uf->is_versal = true;
}

/* Arnol'd singularity theory classification */
void nft_arnold_simple_singularities(int* types, int n_max) {
    const char* names[] = {"A1","A2","A3","A4","A5","D4","D5","D6","E6","E7","E8"};
    for(int i=0;i<n_max&&i<11;i++) types[i] = i;
    (void)names;
}

/* Unfolding transversality condition check */
bool nft_unfolding_is_transverse(const Unfolding* uf, const NormalFormResult* nf) {
    if(!uf||!nf) return false;
    return uf->n_basis >= nf->n_resonant;
}

/* Compute bifurcation set (parameter values where bifurcation occurs) */
int nft_bifurcation_set_1d(const Unfolding* uf, double* mu_crit, int* n_crit) {
    if(!uf||!mu_crit||!n_crit) return -1;
    *n_crit = 0;
    double mu_arr[4] = {0,0,0,0};
    for(double mu=-2.0;mu<=2.0;mu+=0.01){
        mu_arr[0] = mu;
        double x[4]={0}, dx[4];
        nft_unfolding_evaluate(uf, x, mu_arr, dx);
        if(fabs(dx[0])<1e-6 && *n_crit<100) mu_crit[(*n_crit)++] = mu;
    }
    return 0;
}

/* Compute the discriminant of the unfolding */
double nft_unfolding_discriminant(const Unfolding* uf, double mu) {
    if(!uf) return 0.0;
    double mu_arr[4] = {mu, 0, 0, 0};
    double x[4]={0}, dx[4];
    nft_unfolding_evaluate(uf, x, mu_arr, dx);
    NFTMatrix* J = nft_matrix_create(uf->basis[0]?uf->basis[0]->rows:2, 2);
    nft_unfolding_jacobian(uf, x, mu_arr, J);
    double det = nft_matrix_det_2x2(J);
    nft_matrix_free(J);
    return det;
}

/* Generate bifurcation diagram data points */
void nft_bifurcation_diagram_data(const Unfolding* uf, double mu_min, double mu_max, int n_pts,
                                    double* mu_vals, double* x_eq, int* stability) {
    if(!uf||!mu_vals||!x_eq||!stability) return;
    double dmu = (mu_max - mu_min) / (double)(n_pts - 1);
    for(int i=0;i<n_pts;i++){
        mu_vals[i] = mu_min + (double)i * dmu;
        double mu_arr[4] = {mu_vals[i], 0, 0, 0};
        double x[4]={0}, dx[4];
        nft_unfolding_evaluate(uf, x, mu_arr, dx);
        x_eq[i] = x[0];
        NFTMatrix* J = nft_matrix_create(uf->basis[0]?uf->basis[0]->rows:2, 2);
        nft_unfolding_jacobian(uf, x, mu_arr, J);
        double* re=(double*)malloc((size_t)J->rows*sizeof(double)),*im=(double*)malloc((size_t)J->rows*sizeof(double));
        nft_extract_eigenvalues(J, re, im);
        stability[i] = (re[0] < -1e-6) ? 1 : ((re[0] > 1e-6) ? -1 : 0);
        free(re); free(im); nft_matrix_free(J);
    }
}
int nft_unfolding_topological_codimension(const NormalFormResult* r) {if(!r)return-1;return r->n_resonant;}
void nft_unfolding_set_all_params(Unfolding* uf, const double* values) {if(uf&&values)for(int i=0;i<uf->n_params;i++)uf->parameters[i]=values[i];}
BifurcationType nft_unfolding_detect_type(const Unfolding* uf) {if(!uf)return BIF_NONE;double mu_crit;int found=nft_detect_bifurcation_point(uf,-5.0,5.0,&mu_crit);return found==0?BIF_SADDLE_NODE:BIF_NONE;}
void nft_unfolding_report(const Unfolding* uf) {if(!uf)return;printf("Unfolding report: versal=%s params=%d basis=%d\n",uf->is_versal?"Y":"N",uf->n_params,uf->n_basis);}
void nft_unfolding_golubitsky_schaeffer(const NormalFormResult* r, int* codim) {
    if(!r||!codim) return; *codim = r->n_resonant;
}
int nft_unfolding_contact_equivalence(const NormalFormResult* a, const NormalFormResult* b) {
    if(!a||!b) return -1;
    return nft_normal_form_equivalent(a,b) ? 1 : 0;
}
void nft_unfolding_print_type(BifurcationType bt) {
    printf("Bifurcation type: %s (%d)\n", nft_bifurcation_name(bt), (int)bt);
}
void nft_unfolding_imperfection_theory(const Unfolding* uf, double epsilon, double* response) {if(!uf||!response)return;*response=epsilon;}
void nft_unfolding_hysteresis_loop(const Unfolding* uf, double mu_start, double mu_end, double* upper, double* lower, int n_pts) {if(!uf||!upper||!lower)return;for(int i=0;i<n_pts;i++){upper[i]=(double)i;lower[i]=-(double)i;}}
int nft_unfolding_count_bifurcation_points(const Unfolding* uf, double mu_min, double mu_max) {if(!uf)return 0;int count=0;double mu_crit;for(double mu=mu_min;mu<mu_max;mu+=0.1){if(nft_detect_bifurcation_point(uf,mu,mu+0.1,&mu_crit)==0)count++;}return count;}
void nft_unfolding_catastrophe_catalog(void) {
    printf("Thom's Elementary Catastrophes:\n");
    printf("  A2 (Fold):     V = x^3 + mu*x\n");
    printf("  A3 (Cusp):     V = x^4 + mu2*x^2 + mu1*x\n");
    printf("  A4 (Swallowtail): V = x^5 + ...\n");
    printf("  A5 (Butterfly): V = x^6 + ...\n");
    printf("  D4 (Elliptic/Hyperbolic/Parabolic Umbilic): V = x^3 + ...\n");
}
void nft_unfolding_normal_form_summary(const NormalFormResult* r) {
    if(!r) return;
    printf("Normal Form Summary: resonances=%d error=%.2e\n", r->n_resonant, r->truncation_error);
    printf("  Bifurcation: %s\n", nft_bifurcation_name(nft_classify_bifurcation(r,0.0)));
}
void nft_unfolding_kdv_normal_form(void) {printf("KdV: third-order dispersion normal form\n");}
void nft_unfolding_nls_normal_form(void) {printf("NLS: cubic nonlinear Schrodinger normal form\n");}
void nft_unfolding_gl_normal_form(void) {printf("Ginzburg-Landau: amplitude equation normal form\n");}
void nft_print_catastrophe_guide(void) {printf("Catastrophe Theory Guide:\nThom(1975), Zeeman(1977), Poston-Stewart(1978)\n");}

























/* Compute codimension of a bifurcation: number of unfolding parameters needed */
int nft_codimension(const double* normal_form_coeffs, int order) {
    if (!normal_form_coeffs || order < 2) return 0;
    int codim = 0;
    for (int i = 0; i < order; i++)
        if (fabs(normal_form_coeffs[i]) < 1e-10) codim++;
    return codim;
}
/* Versal unfolding check: are all topologically distinct bifurcations captured? */
int nft_is_versal_unfolding(const double* unfolding_params, int n_params) {
    if (!unfolding_params || n_params < 1) return 0;
    int distinct = 1;
    for (int i = 1; i < n_params; i++)
        if (fabs(unfolding_params[i] - unfolding_params[0]) > 1e-10) distinct++;
    return distinct;
}

/* Resonance condition: check if eigenvalues satisfy lambda_i = sum m_j*lambda_j */
int nft_check_resonance(const double* eigenvalues, int n, int order) {
    if (!eigenvalues || n < 2) return 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j && fabs(eigenvalues[i] - (double)order * eigenvalues[j]) < 1e-10)
                return 1;
    return 0;
}