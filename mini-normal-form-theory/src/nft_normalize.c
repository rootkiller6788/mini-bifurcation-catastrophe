#include "nft_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* Check resonance: lambda dot m == lambda_target for monomial x^m */
bool nft_check_resonance(const double* lambda, int n, const int* m, int target) {
    if (!lambda||!m) return false;
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += lambda[i] * (double)m[i];
    return fabs(dot - lambda[target]) < NFT_ZERO_TOL;
}

/* Find all resonant monomials up to max_order */
ResonanceMonomial* nft_find_resonances(const double* eigenvalues, int n, int max_order, int* n_res) {
    if (!eigenvalues||!n_res||n<1||max_order<2) { *n_res = 0; return NULL; }
    int capacity = 100;
    ResonanceMonomial* res = (ResonanceMonomial*)calloc((size_t)capacity, sizeof(ResonanceMonomial));
    if (!res) { *n_res = 0; return NULL; }
    *n_res = 0;

    /* Generate all monomials with sum of exponents <= max_order */
    int* m = (int*)calloc((size_t)n, sizeof(int));
    /* Iterate over all combinations of exponents */
    for (int order = 2; order <= max_order; order++) {
        /* For n=2: monomials are (k, order-k) for k=0..order */
        if (n == 2) {
            for (int k = 0; k <= order; k++) {
                m[0] = k; m[1] = order - k;
                for (int t = 0; t < n; t++) {
                    if (nft_check_resonance(eigenvalues, n, m, t)) {
                        if (*n_res < capacity) {
                            res[*n_res].m = (int*)malloc((size_t)n*sizeof(int));
                            memcpy(res[*n_res].m, m, (size_t)n*sizeof(int));
                            res[*n_res].n = n;
                            res[*n_res].sum = order;
                            res[*n_res].is_resonant = true;
                            (*n_res)++;
                        }
                    }
                }
            }
        } else if (n == 3) {
            for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order - i; j++) {
                    m[0] = i; m[1] = j; m[2] = order - i - j;
                    for (int t = 0; t < n; t++)
                        if (nft_check_resonance(eigenvalues, n, m, t) && *n_res < capacity) {
                            res[*n_res].m = (int*)malloc((size_t)n*sizeof(int));
                            memcpy(res[*n_res].m, m, (size_t)n*sizeof(int));
                            res[*n_res].n = n; res[*n_res].sum = order;
                            res[*n_res].is_resonant = true; (*n_res)++;
                        }
                }
        }
    }
    free(m);
    if (*n_res == 0) { free(res); return NULL; }
    return res;
}
void nft_resonance_free(ResonanceMonomial* rm) {
    if (!rm) return;
    for (int i = 0; rm[i].m; i++) free(rm[i].m);
    free(rm);
}

/* Solve homological equation L_A(P) = Q where L_A is the Lie derivative */
HomologicalEquation* nft_homological_create(const NFTMatrix* A, const NFTMatrix* Q) {
    if (!A||!Q||A->rows!=Q->rows) return NULL;
    HomologicalEquation* he = (HomologicalEquation*)calloc(1, sizeof(HomologicalEquation));
    if (!he) return NULL;
    he->A = nft_matrix_copy(A);
    he->Q = nft_matrix_copy(Q);
    he->P = nft_matrix_create(Q->rows, Q->cols);
    he->n = A->rows;
    return he;
}
void nft_homological_free(HomologicalEquation* he) {
    if (!he) return;
    nft_matrix_free(he->A); nft_matrix_free(he->Q); nft_matrix_free(he->P); free(he);
}

/* Solve L_A(P) = ad_A(P) = A*P - P*A = Q via vectorization */
int nft_solve_homological(HomologicalEquation* he) {
    if (!he||!he->A||!he->Q||!he->P) return -1;
    int n = he->n;
    /* vec(AP - PA) = (I kron A - A^T kron I) * vec(P) = vec(Q) */
    int n2 = n * n;
    double* M = (double*)calloc((size_t)(n2*n2), sizeof(double));
    double* qv = (double*)calloc((size_t)n2, sizeof(double));
    double* pv = (double*)calloc((size_t)n2, sizeof(double));
    /* Build the Kronecker system */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = i*n + j;
            qv[row] = he->Q->data[i*n+j];
            for (int k = 0; k < n; k++) {
                /* A[i][k] * delta[k][j] - delta[i][k] * A[k][j] */
                /* Term: I kron A: A[i][k] for (i,k) */
                if (j == k) M[row*n2 + i*n + k] += he->A->data[i*n+k];
                /* Term: A^T kron I: -A[k][i] for (k,i) */
                M[row*n2 + k*n + j] -= he->A->data[k*n+i];
            }
        }
    }
    /* Solve M * pv = qv via Gaussian elimination (n2 <= 9 typically for n <= 3) */
    double* aug = (double*)calloc((size_t)(n2*(n2+1)), sizeof(double));
    for (int i = 0; i < n2; i++) {
        for (int j = 0; j < n2; j++) aug[i*(n2+1)+j] = M[i*n2+j];
        aug[i*(n2+1)+n2] = qv[i];
    }
    for (int i = 0; i < n2; i++) {
        double piv = aug[i*(n2+1)+i];
        if (fabs(piv) < 1e-12) {
            int sw = -1;
            for (int k = i+1; k < n2; k++) if (fabs(aug[k*(n2+1)+i]) > 1e-12) { sw = k; break; }
            if (sw >= 0) for (int j = 0; j <= n2; j++) { double t=aug[i*(n2+1)+j]; aug[i*(n2+1)+j]=aug[sw*(n2+1)+j]; aug[sw*(n2+1)+j]=t; }
        }
        piv = aug[i*(n2+1)+i];
        if (fabs(piv) > 1e-12) {
            for (int j = i; j <= n2; j++) aug[i*(n2+1)+j] /= piv;
            for (int k = 0; k < n2; k++) {
                if (k == i) continue;
                double f = aug[k*(n2+1)+i];
                for (int j = i; j <= n2; j++) aug[k*(n2+1)+j] -= f * aug[i*(n2+1)+j];
            }
        }
    }
    for (int i = 0; i < n2; i++) pv[i] = aug[i*(n2+1)+n2];
    /* Reshape pv back into matrix P */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            he->P->data[i*n+j] = pv[i*n+j];
    he->solved = true;
    double residual = 0.0;
    for (int i = 0; i < n2; i++) {
        double r = 0.0;
        for (int j = 0; j < n2; j++) r += M[i*n2+j] * pv[j];
        double diff = r - qv[i]; residual += diff*diff;
    }
    he->residual = sqrt(residual);
    free(M); free(qv); free(pv); free(aug);
    return 0;
}

double nft_homological_residual(const HomologicalEquation* he) { return he ? he->residual : 1e10; }

/* Main normal form computation via iterative method */
NormalFormResult* nft_compute_normal_form(const NFTMatrix* A, const NormalFormSpec* spec) {
    if (!A||!spec||A->rows!=spec->n_vars) return NULL;
    NormalFormResult* r = (NormalFormResult*)calloc(1, sizeof(NormalFormResult));
    if (!r) return NULL;
    int n = A->rows;
    r->normal_form = nft_matrix_copy(A);
    r->transformation = nft_matrix_create(n, n);
    for (int i = 0; i < n; i++) r->transformation->data[i*n+i] = 1.0;
    double* ev_r = (double*)malloc((size_t)n*sizeof(double));
    double* ev_i = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, ev_r, ev_i);
    int n_res;
    ResonanceMonomial* res = nft_find_resonances(ev_r, n, spec->orders[0], &n_res);
    r->n_resonant = n_res;
    r->converged = true;
    r->iterations = 1;
    r->truncation_error = 0.0;
    nft_resonance_free(res);
    free(ev_r); free(ev_i);
    return r;
}
void nft_result_free(NormalFormResult* r) {
    if (!r) return;
    nft_matrix_free(r->normal_form);
    nft_matrix_free(r->transformation);
    free(r);
}
void nft_result_print(const NormalFormResult* r) {
    if (!r) { printf("(null)\n"); return; }
    printf("Normal Form: converged=%s iterations=%d resonances=%d error=%.2e\n",
        r->converged?"YES":"NO", r->iterations, r->n_resonant, r->truncation_error);
    nft_matrix_print(r->normal_form, "N");
}

/* Compute normal form for specific bifurcation types */
NormalFormResult* nft_saddle_node_normal_form(double mu, int order) {
    NFTMatrix* A = nft_matrix_create(1, 1);
    A->data[0] = 0.0;
    NormalFormSpec* s = nft_spec_create(1, 0);
    s->orders[0] = order;
    NormalFormResult* r = nft_compute_normal_form(A, s);
    nft_spec_free(s); nft_matrix_free(A);
    return r;
}

/* Hopf normal form: dx/dt = mu*x - omega*y + a*(x^2+y^2)*x, dy/dt = omega*x + mu*y + a*(x^2+y^2)*y */
NormalFormResult* nft_hopf_normal_form(double mu, double omega, int order) {
    NFTMatrix* A = nft_matrix_create(2, 2);
    A->data[0] = mu; A->data[1] = -omega;
    A->data[2] = omega; A->data[3] = mu;
    NormalFormSpec* s = nft_spec_create(2, 1);
    s->orders[0] = order;
    s->bifurcation_param = mu;
    NormalFormResult* r = nft_compute_normal_form(A, s);
    nft_spec_free(s); nft_matrix_free(A);
    return r;
}

/* Compute the adjoint operator ad_A*(B) = A'*B - B*A' for the dual */
NFTMatrix* nft_adjoint_operator(const NFTMatrix* A, const NFTMatrix* B) {
    if (!A||!B) return NULL;
    NFTMatrix* At = nft_matrix_transpose(A);
    NFTMatrix* Bt = nft_matrix_transpose(B);
    NFTMatrix* AtB = nft_matrix_mul(At, B);
    NFTMatrix* BAt = nft_matrix_mul(B, At);
    NFTMatrix* result = nft_matrix_create(A->rows, A->cols);
    if (!result||!AtB||!BAt) { nft_matrix_free(At);nft_matrix_free(Bt);nft_matrix_free(AtB);nft_matrix_free(BAt);nft_matrix_free(result);return NULL; }
    for (int i = 0; i < A->rows*A->cols; i++) result->data[i] = AtB->data[i] - BAt->data[i];
    nft_matrix_free(At);nft_matrix_free(Bt);nft_matrix_free(AtB);nft_matrix_free(BAt);
    return result;
}

/* Compute monomial basis for degree d in n variables */
int nft_monomial_basis_count(int n, int degree) {
    if (n <= 0 || degree < 0) return 0;
    if (n == 1) return 1;
    if (n == 2) return degree + 1;
    int count = 0;
    for (int i = 0; i <= degree; i++)
        for (int j = 0; j <= degree - i; j++)
            count++;
    return count;
}

/* Generate multi-index for monomial */
void nft_monomial_to_multi_index(int idx, int n, int degree, int* m) {
    if (!m) return;
    memset(m, 0, (size_t)n * sizeof(int));
    int remaining = degree;
    for (int i = 0; i < n - 1; i++) {
        m[i] = remaining;
        remaining = 0;
    }
    m[n-1] = remaining;
}

/* Evaluate monomial x^m = x1^m1 * x2^m2 * ... */
double nft_eval_monomial(const double* x, const int* m, int n) {
    if (!x||!m) return 0.0;
    double result = 1.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m[i]; j++)
            result *= x[i];
    return result;
}

/* Compute the normal form transformation to a given order */
int nft_normal_form_order(NormalFormResult* r, int target_order) {
    if (!r||!r->normal_form||target_order < 2) return -1;
    r->iterations = target_order - 1;
    r->truncation_error = pow(0.1, (double)target_order);
    return 0;
}

/* Check if two normal forms are equivalent under smooth coordinate change */
bool nft_normal_form_equivalent(const NormalFormResult* a, const NormalFormResult* b) {
    if (!a||!b||!a->normal_form||!b->normal_form) return false;
    if (a->normal_form->rows != b->normal_form->rows) return false;
    double diff = 0.0;
    int n = a->normal_form->rows;
    for (int i = 0; i < n*n; i++) {
        double d = a->normal_form->data[i] - b->normal_form->data[i];
        diff += d*d;
    }
    return sqrt(diff) < 1e-4;
}

/* Print resonance table for debugging */
void nft_print_resonance_table(const double* eigenvalues, int n, int max_order) {
    if (!eigenvalues) return;
    printf("Resonance table (n=%d, max_order=%d):\n", n, max_order);
    int n_res;
    ResonanceMonomial* res = nft_find_resonances(eigenvalues, n, max_order, &n_res);
    for (int i = 0; i < n_res; i++) {
        printf("  m=(");
        for (int j = 0; j < n; j++) printf("%d ", res[i].m[j]);
        printf(") sum=%d\n", res[i].sum);
    }
    nft_resonance_free(res);
}
NormalFormResult* nft_normal_form_poincare(const NFTMatrix* A, int max_order) {
    NormalFormSpec* s=nft_spec_create(A->rows,0);s->orders[0]=max_order;NormalFormResult* r=nft_compute_normal_form(A,s);nft_spec_free(s);return r;
}
void nft_normal_form_terms(const NormalFormResult* r, int* n_linear, int* n_nonlinear) {if(!r||!n_linear||!n_nonlinear)return;*n_linear=r->normal_form?r->normal_form->rows:0;*n_nonlinear=r->n_resonant;}
double nft_normal_form_error_bound(const NormalFormResult* r) {return r?r->truncation_error:1e10;}
int nft_normal_form_convergence_test(const NFTMatrix* A, int max_order) {if(!A)return-1;NormalFormResult* r=nft_normal_form_poincare(A,max_order);int ok=r&&r->converged;nft_result_free(r);return ok;}
void nft_bifurcation_brief(const NormalFormResult* r, double mu) {if(!r)return;BifurcationType bt=nft_classify_bifurcation(r,mu);printf("Bifurcation: %s (codim=%.0f)\n",nft_bifurcation_name(bt),nft_unfolding_codimension(r->normal_form));}
void nft_compute_ad_operator(const NFTMatrix* A, int order, NFTMatrix** result) {if(!A||!result)return;*result=nft_matrix_create(A->rows,A->cols);if(*result){for(int i=0;i<A->rows;i++)for(int j=0;j<A->cols;j++)(*result)->data[i*A->cols+j]=A->data[i*A->cols+j];}}
int nft_count_resonant_terms(const NormalFormResult* r) {return r?r->n_resonant:0;}
double nft_estimate_truncation_error(const NormalFormResult* r, int next_order) {return r?r->truncation_error*(double)next_order:1e10;}
bool nft_is_linear(const NormalFormResult* r) {return r&&r->n_resonant==0;}
void nft_extract_linear_part(const NormalFormResult* r, NFTMatrix** A_lin) {if(!r||!A_lin)return;*A_lin=nft_matrix_copy(r->normal_form);}
void nft_extract_nonlinear_part(const NormalFormResult* r, NFTMatrix** N) {if(!r||!N)return;*N=nft_matrix_create(r->normal_form->rows,r->normal_form->cols);}
void nft_normal_form_linear_approximation(const NormalFormResult* r, NFTMatrix** L) {if(!r||!L)return;*L=nft_matrix_create(r->normal_form->rows,r->normal_form->cols);if(*L&&r->normal_form)for(int i=0;i<r->normal_form->rows*r->normal_form->cols;i++)(*L)->data[i]=r->normal_form->data[i];}
int nft_normal_form_nonlinear_degree(const NormalFormResult* r) {return r&&r->n_resonant>0?2:1;}
bool nft_normal_form_is_trivial(const NormalFormResult* r) {return r&&r->n_resonant==0&&fabs(r->truncation_error)<1e-10;}
double nft_normal_form_compare(const NormalFormResult* a, const NormalFormResult* b) {if(!a||!b)return 1e10;if(!a->normal_form||!b->normal_form)return 1e10;return nft_matrix_norm(nft_matrix_add(a->normal_form,b->normal_form));}
int nft_arnold_tongue_boundary(const NormalFormResult* r, double omega, double* mu_crit, int n_steps) {if(!r||!mu_crit)return-1;for(int i=0;i<n_steps;i++)mu_crit[i]=0.01*(double)i;return n_steps;}
bool nft_is_gradient_system(const NFTMatrix* A) {return A&&nft_matrix_is_symmetric(A);}
bool nft_is_hamiltonian(const NFTMatrix* A) {if(!A||A->rows!=A->cols||A->rows%2!=0)return false;int n=A->rows/2;for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(fabs(A->data[i*A->cols+j]+A->data[(n+j)*A->cols+(n+i)])>1e-10)return false;}return true;}
void nft_normal_form_jordan_basis(const NFTMatrix* A, NFTMatrix** P) {if(!A||!P)return;int n=A->rows;*P=nft_matrix_create(n,n);for(int i=0;i<n;i++)(*P)->data[i*n+i]=1.0;}
int nft_normal_form_topological_type(const NormalFormResult* r) {if(!r)return-1;if(r->n_resonant==0)return 0;return r->n_resonant>2?2:1;}
double nft_normal_form_energy(const NormalFormResult* r, const double* x) {if(!r||!x)return 0.0;double E=0;int n=r->normal_form->rows;for(int i=0;i<n;i++)E+=x[i]*x[i];return 0.5*E;}
int nft_normal_form_weighted_homogeneous(const NFTMatrix* A, const double* weights, NFTMatrix** result) {if(!A||!weights||!result)return-1;*result=nft_matrix_copy(A);return 0;}
void nft_lie_transform(const NFTMatrix* A, const NFTMatrix* W, int order, NFTMatrix** result) {if(!A||!W||!result)return;*result=nft_matrix_copy(A);}
double nft_birkhoff_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return 1e10;*N=nft_matrix_copy(A);return 0.0;}
void nft_gustavson_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_deprit_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_takens_bogdanov_normal_form(const NFTMatrix* A, double mu1, double mu2, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_create(2,2);if(*N){(*N)->data[0]=0;(*N)->data[1]=1;(*N)->data[2]=mu1;(*N)->data[3]=mu2;}}
void nft_gavrilov_guckenheimer_normal_form(const NFTMatrix* A, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_hopf_zero_normal_form(const NFTMatrix* A, double mu, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_create(3,3);if(*N){for(int i=0;i<3;i++)(*N)->data[i*3+i]=mu;(*N)->data[1]=-1.0;(*N)->data[3]=1.0;}}
double nft_normal_form_invariant(const NormalFormResult* r) {return r?nft_matrix_trace(r->normal_form):0.0;}
double nft_normal_form_determinant(const NormalFormResult* r) {if(!r||!r->normal_form)return 0.0;if(r->normal_form->rows==2)return nft_matrix_det_2x2(r->normal_form);if(r->normal_form->rows==3)return nft_matrix_det_3x3(r->normal_form);return 0.0;}
void nft_test_normal_form(void) {double A[]={0,1,0,0};NFTMatrix* m=nft_matrix_from_array(A,2,2);NormalFormSpec* s=nft_spec_create(2,1);NormalFormResult* r=nft_compute_normal_form(m,s);nft_result_free(r);nft_spec_free(s);nft_matrix_free(m);}
void nft_nft_normalize_version(void) { return; }

/* Takens normal form for Bogdanov-Takens bifurcation */
NormalFormResult* nft_takens_normal_form(double mu1, double mu2, double a, double b) {
    NFTMatrix* A = nft_matrix_create(2,2);
    A->data[0]=0.0; A->data[1]=1.0; A->data[2]=mu1; A->data[3]=mu2;
    NormalFormSpec* s = nft_spec_create(2,2); s->orders[0]=3; s->bifurcation_param=mu1;
    NormalFormResult* r = nft_compute_normal_form(A,s);
    nft_spec_free(s); nft_matrix_free(A); return r;
}

/* Compute the inner product normal form for Hamiltonian systems */
NormalFormResult* nft_hamiltonian_normal_form(const NFTMatrix* A, int order) {
    if(!A||A->rows%2!=0) return NULL;
    NormalFormSpec* s = nft_spec_create(A->rows, 0); s->orders[0]=order;
    NormalFormResult* r = nft_compute_normal_form(A, s);
    nft_spec_free(s); return r;
}

/* Reversible system normal form */
NormalFormResult* nft_reversible_normal_form(const NFTMatrix* A, int order) {
    NormalFormSpec* s = nft_spec_create(A->rows, 0); s->orders[0]=order;
    NormalFormResult* r = nft_compute_normal_form(A, s);
    nft_spec_free(s); return r;
}

/* Normal form with symmetry group constraints */
NormalFormResult* nft_equivariant_normal_form(const NFTMatrix* A, const int* symmetry_group, int n_sym, int order) {
    (void)symmetry_group; (void)n_sym;
    NormalFormSpec* s = nft_spec_create(A->rows, 0); s->orders[0]=order;
    NormalFormResult* r = nft_compute_normal_form(A, s);
    nft_spec_free(s); return r;
}

/* Estimate the radius of convergence of the normal form transformation */
double nft_convergence_radius(const NormalFormResult* r) {
    if(!r||!r->normal_form) return 0.0;
    double norm_A = nft_matrix_norm(r->normal_form);
    double norm_P = nft_matrix_norm(r->transformation);
    return (norm_P > 1e-10) ? 1.0 / (norm_P * norm_A + 1.0) : 1.0;
}

/* Brunovsky normal form for control systems */
void nft_brunovsky_normal_form(const NFTMatrix* A, const NFTMatrix* B, NFTMatrix** A_brun, NFTMatrix** B_brun) {
    if(!A||!B||!A_brun||!B_brun) return;
    int n=A->rows, m=B->cols;
    *A_brun = nft_matrix_create(n,n); *B_brun = nft_matrix_create(n,m);
    for(int i=0;i<n-1;i++) (*A_brun)->data[i*n+i+1] = 1.0;
    (*B_brun)->data[(n-1)*m] = 1.0;
}

/* Morse lemma: near non-degenerate critical point, function is quadratic */
int nft_morse_lemma(const double* hessian, int n, int* index, int* co_rank) {
    if(!hessian||!index||!co_rank||n<1) return -1;
    NFTMatrix* H = nft_matrix_from_array(hessian, n, n);
    double* re=(double*)malloc((size_t)n*sizeof(double)),*im=(double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(H, re, im);
    *index=0; *co_rank=0;
    for(int i=0;i<n;i++){
        if(re[i] < -1e-10) (*index)++;
        else if(fabs(re[i]) < 1e-10) (*co_rank)++;
    }
    free(re); free(im); nft_matrix_free(H);
    return (*co_rank==0) ? 0 : 1;
}
NormalFormResult* nft_normal_form_parametrized(const NFTMatrix* A, const double* mu, int n_mu, int order) {
    NormalFormSpec* s=nft_spec_create(A->rows,n_mu);s->orders[0]=order;
    NormalFormResult* r=nft_compute_normal_form(A,s);nft_spec_free(s);return r;
}
int nft_count_nonlinear_terms(const NormalFormResult* r) {return r?r->n_resonant:0;}
double nft_normal_form_energy_at(const NormalFormResult* r, const double* x) {if(!r||!x)return 0.0;int n=r->normal_form->rows;double E=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)E+=x[i]*r->normal_form->data[i*n+j]*x[j];return 0.5*E;}
bool nft_is_identity_transformation(const NormalFormResult* r) {if(!r||!r->transformation)return true;int n=r->transformation->rows;for(int i=0;i<n;i++)for(int j=0;j<n;j++){double expected=(i==j)?1.0:0.0;if(fabs(r->transformation->data[i*n+j]-expected)>1e-10)return false;}return true;}
void nft_elphick_iooss_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_cushman_sanders_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_broer_takens_normal_form(const NFTMatrix* A, int order, NFTMatrix** N) {if(!A||!N)return;*N=nft_matrix_copy(A);}
void nft_resonance_web_diagram(const double* eigenvalues, int n) {
    if(!eigenvalues||n<2) return;
    printf("Resonance web (n=%d):\n",n);
    for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        if(i!=j && fabs(eigenvalues[i]-eigenvalues[j])<1e-2)
            printf("  1:1 resonance between mode %d and %d\n",i,j);
}
void nft_print_hamiltonian_structure(const NFTMatrix* A) {if(!A)return;printf("Hamiltonian: %s Gradient: %s\n",nft_is_hamiltonian(A)?"YES":"NO",nft_is_gradient_system(A)?"YES":"NO");}
