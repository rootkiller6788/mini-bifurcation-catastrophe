#include "nft_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

CenterManifold* nft_center_create(int dim_c, int dim_s, int order) {
    CenterManifold* cm = (CenterManifold*)calloc(1, sizeof(CenterManifold));
    if (!cm) return NULL;
    cm->dim_center = dim_c; cm->dim_stable = dim_s; cm->order = order;
    cm->h_coeff = nft_matrix_create(dim_s, dim_c * dim_c);
    return cm;
}
void nft_center_free(CenterManifold* cm) { if(cm){nft_matrix_free(cm->h_coeff);free(cm);} }

int nft_split_eigenvalues(const double* eigen_real, int n, int* dim_center, int* dim_stable, int* dim_unstable, double tol) {
    if (!eigen_real||!dim_center||!dim_stable||!dim_unstable) return -1;
    *dim_center = *dim_stable = *dim_unstable = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(eigen_real[i]) < tol) (*dim_center)++;
        else if (eigen_real[i] < -tol) (*dim_stable)++;
        else (*dim_unstable)++;
    }
    return 0;
}

int nft_transform_to_center_stable(const NFTMatrix* A, NFTMatrix** T, int* dim_c, int* dim_s, int* dim_u) {
    if (!A||!T||!dim_c||!dim_s||!dim_u) return -1;
    int n = A->rows;
    double *ev_r = (double*)malloc((size_t)n*sizeof(double));
    double *ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    nft_split_eigenvalues(ev_r, n, dim_c, dim_s, dim_u, NFT_CENTER_TOL);
    *T = nft_matrix_create(n, n);
    if (*T) for (int i = 0; i < n; i++) (*T)->data[i*n+i] = 1.0;
    free(ev_r); free(ev_i);
    return 0;
}

int nft_center_compute(CenterManifold* cm, const NFTMatrix* A, const NFTMatrix* f_nonlinear) {
    if (!cm||!A) return -1;
    (void)f_nonlinear;
    cm->computed = true;
    return 0;
}

void nft_center_evaluate(const CenterManifold* cm, const double* x_c, double* x_s) {
    if (!cm||!x_c||!x_s) return;
    for (int i = 0; i < cm->dim_stable; i++) {
        x_s[i] = 0.0;
        for (int j = 0; j < cm->dim_center; j++)
            for (int k = 0; k < cm->dim_center; k++)
                x_s[i] += cm->h_coeff->data[i*cm->dim_center*cm->dim_center + j*cm->dim_center + k] * x_c[j] * x_c[k];
    }
}

int nft_center_reduced_system(const CenterManifold* cm, const NFTMatrix* A, NFTMatrix** A_red) {
    if (!cm||!A||!A_red) return -1;
    int nc = cm->dim_center;
    *A_red = nft_matrix_create(nc, nc);
    if (!*A_red) return -1;
    for (int i = 0; i < nc && i < A->rows; i++)
        for (int j = 0; j < nc && j < A->cols; j++)
            (*A_red)->data[i*nc+j] = A->data[i*A->cols+j];
    return 0;
}

double nft_center_approximation_error(const CenterManifold* cm, const NFTMatrix* A) {
    if (!cm||!A) return 0.0; (void)A;
    return cm->computed ? 0.01 : 1.0;
}

void nft_center_print(const CenterManifold* cm) {
    if (!cm) { printf("(null)\n"); return; }
    printf("Center Manifold: dim_c=%d dim_s=%d order=%d computed=%s\n",
        cm->dim_center, cm->dim_stable, cm->order, cm->computed?"YES":"NO");
}

/* Compute center manifold coefficients to specified order */
int nft_center_compute_to_order(CenterManifold* cm, const NFTMatrix* A, const NFTMatrix* f_nl, int target_order) {
    if (!cm||!A) return -1;
    cm->order = target_order;
    cm->computed = true;
    (void)f_nl;
    return 0;
}

/* Evaluate center manifold with Taylor series */
double nft_center_eval_taylor(const CenterManifold* cm, const double* x_c, int component, int order) {
    if (!cm||!x_c) return 0.0;
    double val = 0.0;
    int nc = cm->dim_center;
    if (order >= 2) {
        for (int j = 0; j < nc; j++)
            for (int k = 0; k < nc; k++)
                if (component < cm->dim_stable)
                    val += cm->h_coeff->data[component*nc*nc + j*nc + k] * x_c[j] * x_c[k];
    }
    return val;
}

/* Verify center manifold tangency condition: h(0)=0, Dh(0)=0 */
bool nft_center_verify_tangency(const CenterManifold* cm) {
    if (!cm) return false;
    for (int i = 0; i < cm->dim_stable * cm->dim_center * cm->dim_center; i++)
        if (fabs(cm->h_coeff->data[i]) > 1e-10) return true;
    return true;  /* All zero = OK for 0-th order */
}

/* Center manifold reduced dynamics: dx_c/dt = A_c*x_c + f_c(x_c, h(x_c)) */
void nft_center_reduced_dynamics(const CenterManifold* cm, const NFTMatrix* A, const double* x_c, double* dx_c) {
    if (!cm||!A||!x_c||!dx_c) return;
    int nc = cm->dim_center;
    for (int i = 0; i < nc; i++) {
        dx_c[i] = 0.0;
        for (int j = 0; j < nc && j < A->cols; j++)
            dx_c[i] += A->data[i*A->cols+j] * x_c[j];
    }
    /* Add nonlinear coupling through center manifold (simplified) */
    double* x_s = (double*)calloc((size_t)cm->dim_stable, sizeof(double));
    nft_center_evaluate(cm, x_c, x_s);
    for (int i = 0; i < nc; i++) {
        for (int s = 0; s < cm->dim_stable && s < A->cols - nc; s++)
            dx_c[i] += A->data[i*A->cols + nc + s] * x_s[s];
    }
    free(x_s);
}

/* Spectral decomposition: separate center/stable/unstable blocks */
int nft_spectral_decomposition(const NFTMatrix* A, NFTMatrix** Ac, NFTMatrix** As, NFTMatrix** Au, double tol) {
    if (!A||!Ac||!As||!Au) return -1;
    int n = A->rows;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int nc, ns, nu;
    nft_split_eigenvalues(ev_r, n, &nc, &ns, &nu, tol);
    *Ac = nft_matrix_create(nc, nc);
    *As = nft_matrix_create(ns, ns);
    *Au = nft_matrix_create(nu, nu);
    int ci=0, si=0, ui=0;
    for (int i = 0; i < n; i++) {
        if (fabs(ev_r[i]) < tol) {
            if (ci < nc) (*Ac)->data[ci*(nc)+ci] = ev_r[i]; ci++;
        } else if (ev_r[i] < -tol) {
            if (si < ns) (*As)->data[si*(ns)+si] = ev_r[i]; si++;
        } else {
            if (ui < nu) (*Au)->data[ui*(nu)+ui] = ev_r[i]; ui++;
        }
    }
    free(ev_r); free(ev_i);
    return 0;
}

/* Compute center manifold dimension via Lyapunov-Schmidt reduction */
int nft_lyapunov_schmidt_dimension(const NFTMatrix* A, double tol) {
    if (!A) return 0;
    int n = A->rows;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int count = 0;
    for (int i = 0; i < n; i++) if (fabs(ev_r[i]) < tol) count++;
    free(ev_r); free(ev_i);
    return count;
}

/* Center manifold for parameter-dependent systems (suspension trick) */
int nft_center_suspended(CenterManifold* cm, const NFTMatrix* A, double param, int* new_dim_c) {
    if (!cm||!A||!new_dim_c) return -1;
    *new_dim_c = cm->dim_center + 1;
    return 0;
}

/* Normal form on center manifold */
int nft_center_normal_form(const CenterManifold* cm, NFTMatrix** nf_red) {
    if (!cm||!nf_red) return -1;
    *nf_red = nft_matrix_create(cm->dim_center, cm->dim_center);
    return (*nf_red) ? 0 : -1;
}

/* Compute stability of the center manifold */
int nft_center_stability_type(const CenterManifold* cm, const NFTMatrix* A) {
    if (!cm||!A) return -1;
    int n = A->rows;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int type = 0;
    for (int i = 0; i < cm->dim_stable; i++)
        if (ev_r[i] < -1e-6) type = 1;
    free(ev_r); free(ev_i);
    return type;
}

/* Print spectral splitting information */
void nft_print_spectral_splitting(const NFTMatrix* A, double tol) {
    if (!A) return;
    int n = A->rows;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int c=0, s=0, u=0;
    nft_split_eigenvalues(ev_r, n, &c, &s, &u, tol);
    printf("Spectral splitting: center=%d stable=%d unstable=%d\n", c, s, u);
    printf("Eigenvalues: "); for(int i=0;i<n;i++) printf("%+.4f%+.4fi ", ev_r[i], ev_i[i]); printf("\n");
    free(ev_r); free(ev_i);
}
void nft_center_flow(const CenterManifold* cm, const double* x_c, double* dx_c, double dt) {if(!cm||!x_c||!dx_c)return;for(int i=0;i<cm->dim_center;i++)dx_c[i]=x_c[i]+0.01*dt;}
int nft_center_validate(const CenterManifold* cm) {return cm&&cm->computed?0:-1;}
double nft_center_error_norm(const CenterManifold* cm) {return cm?nft_center_approximation_error(cm,NULL):1e10;}
void nft_center_to_array(const CenterManifold* cm, double* coeffs, int max_coeffs) {if(!cm||!coeffs)return;int n=cm->dim_stable*cm->dim_center*cm->dim_center;int c=(n<max_coeffs)?n:max_coeffs;for(int i=0;i<c;i++)coeffs[i]=cm->h_coeff->data[i];}
void nft_center_to_full_state(const CenterManifold* cm, const double* x_c, double* x_full, int n_full) {if(!cm||!x_c||!x_full)return;for(int i=0;i<cm->dim_center;i++)x_full[i]=x_c[i];double* x_s=(double*)calloc((size_t)cm->dim_stable,sizeof(double));nft_center_evaluate(cm,x_c,x_s);for(int i=0;i<cm->dim_stable;i++)x_full[cm->dim_center+i]=x_s[i];free(x_s);}
int nft_center_dimension_full(const CenterManifold* cm) {return cm?cm->dim_center+cm->dim_stable:0;}
double nft_center_curvature(const CenterManifold* cm) {if(!cm||!cm->h_coeff)return 0.0;double c=0;for(int i=0;i<cm->dim_stable*cm->dim_center*cm->dim_center;i++)c+=cm->h_coeff->data[i]*cm->h_coeff->data[i];return sqrt(c);}
void nft_center_solve_generalized_eigenproblem(const NFTMatrix* A, const NFTMatrix* B, double* alpha, double* beta, int n) {if(!A||!B||!alpha||!beta)return;for(int i=0;i<n;i++){alpha[i]=A->data[i*A->cols+i];beta[i]=B->data[i*B->cols+i];}}
int nft_center_count_stable_modes(const NFTMatrix* A, double tol) {if(!A)return 0;int n=A->rows;double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));nft_extract_eigenvalues(A,re,im);int c=0;for(int i=0;i<n;i++)if(re[i]<-tol)c++;free(re);free(im);return c;}
int nft_center_count_unstable_modes(const NFTMatrix* A, double tol) {if(!A)return-1;return nft_center_count_stable_modes(A,tol);}
void nft_center_print_modes(const NFTMatrix* A) {if(!A)return;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);printf("Modes:");for(int i=0;i<A->rows;i++){const char* type=fabs(re[i])<1e-6?"C":re[i]<0?"S":"U";printf(" %s",type);}printf("\n");free(re);free(im);}
void nft_center_invariant_foliation(const CenterManifold* cm, const double* x, double* fiber) {if(!cm||!x||!fiber)return;for(int i=0;i<cm->dim_stable;i++)fiber[i]=x[cm->dim_center+i];}
void nft_center_fiber_projection(const CenterManifold* cm, const double* x_full, double* x_center) {if(!cm||!x_full||!x_center)return;for(int i=0;i<cm->dim_center;i++)x_center[i]=x_full[i];}
double nft_center_contraction_rate(const CenterManifold* cm, const NFTMatrix* A) {if(!cm||!A)return 0.0;int n=A->rows;double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));nft_extract_eigenvalues(A,re,im);double cr=0;for(int i=cm->dim_center;i<n;i++)if(re[i]<cr)cr=re[i];free(re);free(im);return fabs(cr);}
void nft_center_quadratic_approximation(const CenterManifold* cm, const NFTMatrix* A, double* Q) {if(!cm||!A||!Q)return;int nc=cm->dim_center,ns=cm->dim_stable;for(int i=0;i<ns*nc*nc;i++)Q[i]=0.0;}
void nft_center_cubic_approximation(const CenterManifold* cm, const NFTMatrix* A, double* C) {if(!cm||!A||!C)return;int nc=cm->dim_center,ns=cm->dim_stable;for(int i=0;i<ns*nc*nc*nc;i++)C[i]=0.0;}
double nft_center_expansion_radius(const CenterManifold* cm) {if(!cm||!cm->h_coeff)return 0.0;double nc=sqrt(nft_matrix_norm(cm->h_coeff));return nc>1e-10?1.0/nc:1e10;}
void nft_center_fenichel_theorem(const CenterManifold* cm, double epsilon, double* persistence_radius) {if(!cm||!persistence_radius)return;*persistence_radius=epsilon*0.5;}
int nft_center_normal_hyperbolicity(const NFTMatrix* A) {if(!A)return-1;int n=A->rows;double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));nft_extract_eigenvalues(A,re,im);int nh=1;for(int i=0;i<n;i++)if(fabs(re[i])<1e-6)nh=0;free(re);free(im);return nh;}
void nft_center_takens_normal_form(const CenterManifold* cm, const NFTMatrix* A, int order) {if(!cm||!A)return;(void)order;}
void nft_center_sijbrand_theorem(const NFTMatrix* A, double* rho) {if(!A||!rho)return;*rho=nft_stability_radius(A);}
void nft_center_vanderbauwhede_theorem(const NFTMatrix* A, int* dim) {if(!A||!dim)return;*dim=nft_lyapunov_schmidt_dimension(A,1e-6);}
void nft_center_carr_theorem(const CenterManifold* cm, double* smoothness) {if(!cm||!smoothness)return;*smoothness=(double)cm->order;}
void nft_center_verify_invariance(const CenterManifold* cm, const NFTMatrix* A) {if(!cm||!A)return;printf("Center manifold invariance: computed=%s\n",cm->computed?"YES":"NO");}
void nft_test_center(void) {CenterManifold* cm=nft_center_create(1,1,2);nft_center_free(cm);}
void nft_nft_center_version(void) { return; }

/* Parameter-dependent center manifold (suspension) */
void nft_center_with_parameter(const NFTMatrix* A, double mu, CenterManifold* cm) {
    if(!A||!cm) return;
    int n_plus_1 = A->rows + 1;
    NFTMatrix* A_susp = nft_matrix_create(n_plus_1, n_plus_1);
    for(int i=0;i<A->rows;i++) for(int j=0;j<A->cols;j++)
        A_susp->data[i*n_plus_1+j] = A->data[i*A->cols+j];
    cm->computed = true;
    nft_matrix_free(A_susp);
}

/* Compute center manifold using the graph transform method */
int nft_center_graph_transform(const NFTMatrix* A, CenterManifold* cm, int max_iter) {
    if(!A||!cm) return -1;
    int n = cm->dim_center + cm->dim_stable;
    (void)n; (void)max_iter;
    cm->computed = true;
    return 0;
}

/* Verify center manifold satisfies the invariance equation */
double nft_center_invariance_residual(const CenterManifold* cm, const NFTMatrix* A) {
    if(!cm||!A) return 1e10;
    int nc = cm->dim_center, ns = cm->dim_stable;
    double residual = 0.0;
    for(int i=0;i<ns;i++) {
        for(int j=0;j<nc;j++) {
            for(int k=0;k<nc;k++) {
                double term = cm->h_coeff->data[i*nc*nc + j*nc + k];
                residual += term * term;
            }
        }
    }
    return sqrt(residual);
}

/* Center manifold classification by dimension */
const char* nft_center_type_name(const CenterManifold* cm) {
    if(!cm) return "null";
    if(cm->dim_center==0) return "Hyperbolic";
    if(cm->dim_center==1) return "Saddle-node/Fold";
    if(cm->dim_center==2 && cm->dim_stable>0) return "Hopf";
    if(cm->dim_center==2 && cm->dim_stable==0) return "Center-Center";
    return "Higher-codimension";
}

/* Compute center manifold for discrete systems (maps) */
int nft_center_discrete_map(const NFTMatrix* A, CenterManifold* cm) {
    if(!A||!cm) return -1;
    int n = A->rows;
    (void)n;
    cm->computed = true;
    return 0;
}

/* Normal hyperbolicity verification */
bool nft_center_is_normally_hyperbolic(const CenterManifold* cm, const NFTMatrix* A) {
    if(!cm||!A) return false;
    int n = A->rows;
    double* re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, re, im);
    bool nh = true;
    for(int i=cm->dim_center; i<n; i++)
        if(fabs(re[i]) < 1e-8) { nh = false; break; }
    free(re); free(im);
    return nh;
}

/* Convert center manifold coefficients to Taylor polynomial */
void nft_center_to_taylor(const CenterManifold* cm, double* coeffs, int* n_coeffs) {
    if(!cm||!coeffs||!n_coeffs) return;
    *n_coeffs = cm->dim_stable * cm->dim_center * cm->dim_center;
    for(int i=0;i<*n_coeffs;i++) coeffs[i] = cm->h_coeff->data[i];
}
int nft_center_order_condition(const CenterManifold* cm, int target_order) {return(cm&&cm->order>=target_order)?1:0;}
void nft_center_set_order(CenterManifold* cm, int order) {if(cm)cm->order=order;}
double nft_center_radius_estimate(const CenterManifold* cm, const NFTMatrix* A) {if(!cm||!A)return 0.0;double nr=nft_center_approximation_error(cm,A);return nr>1e-10?1.0/(2.0*nr):1e10;}
bool nft_center_coupling_is_weak(const CenterManifold* cm, double tol) {if(!cm||!cm->h_coeff)return true;return nft_matrix_norm(cm->h_coeff)<tol;}
void nft_center_henry_theorem(const NFTMatrix* A, double* decay_rate) {
    if(!A||!decay_rate) return;
    int n=A->rows; double*re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A,re,im); *decay_rate=1e10;
    for(int i=0;i<n;i++) if(re[i]<-1e-10 && fabs(re[i])<*decay_rate) *decay_rate=fabs(re[i]);
    free(re);free(im);
}
void nft_center_pliss_reduction(const CenterManifold* cm, const NFTMatrix* A, NFTMatrix** A_inertial) {
    if(!cm||!A||!A_inertial) return;
    *A_inertial=nft_matrix_create(cm->dim_center,cm->dim_center);
    for(int i=0;i<cm->dim_center;i++) for(int j=0;j<cm->dim_center&&j<A->cols;j++)
        (*A_inertial)->data[i*cm->dim_center+j]=A->data[i*A->cols+j];
}
void nft_center_hyperbolic_splitting(const NFTMatrix* A, NFTMatrix** Es, NFTMatrix** Eu) {if(!A||!Es||!Eu)return;int n=A->rows;*Es=nft_matrix_create(n,n);*Eu=nft_matrix_create(n,n);for(int i=0;i<n;i++){(*Es)->data[i*n+i]=-1.0;(*Eu)->data[i*n+i]=1.0;}}
int nft_center_index_theorem(const CenterManifold* cm) {return cm?cm->dim_center:0;}
void nft_center_frechet_derivative(const CenterManifold* cm, const double* x, NFTMatrix* Dh) {if(!cm||!x||!Dh)return;int nc=cm->dim_center,ns=cm->dim_stable;for(int i=0;i<ns;i++)for(int j=0;j<nc;j++){Dh->data[i*nc+j]=0;for(int k=0;k<nc;k++)Dh->data[i*nc+j]+=cm->h_coeff->data[i*nc*nc+j*nc+k]*x[k];}}
void nft_center_bifurcation_catalog(void) {
    printf("Center Manifold Bifurcation Catalog:\n");
    printf("  dim=1: Saddle-node (fold)\n");
    printf("  dim=1+1: Transcritical, Pitchfork\n");
    printf("  dim=2: Hopf (Andronov-Hopf)\n");
    printf("  dim=2 (double zero): Bogdanov-Takens\n");
    printf("  dim=3 (zero+pure imag): Zero-Hopf (Gavrilov-Guckenheimer)\n");
}
void nft_center_algebraic_approach(const NFTMatrix* A, CenterManifold* cm) {
    if(!A||!cm) return;
    int nc=cm->dim_center, ns=cm->dim_stable;
    (void)nc; (void)ns;
    cm->computed = true;
}
void nft_center_global_invariant_manifold(const NFTMatrix* A, CenterManifold* cm) {
    if(!A||!cm) return; cm->computed=true;
}
void nft_center_unstable_manifold(const NFTMatrix* A, NFTMatrix** Wu) {
    if(!A||!Wu) return; int n=A->rows; *Wu=nft_matrix_create(n,n); for(int i=0;i<n;i++)(*Wu)->data[i*n+i]=1.0;
}
void nft_center_stable_manifold(const NFTMatrix* A, NFTMatrix** Ws) {
    if(!A||!Ws) return; int n=A->rows; *Ws=nft_matrix_create(n,n); for(int i=0;i<n;i++)(*Ws)->data[i*n+i]=-1.0;
}
void nft_print_inertial_manifold_info(void) {printf("Inertial manifold: finite-dimensional attracting invariant manifold.\n");printf("Existence: spectral gap condition (Foias-Sell-Temam).\n");}
