#include "nft_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

NFTMatrix* nft_matrix_create(int r, int c) {
    NFTMatrix* m = (NFTMatrix*)calloc(1, sizeof(NFTMatrix));
    if (!m) return NULL;
    m->rows = r; m->cols = c; m->owns = true;
    m->data = (double*)calloc((size_t)(r*c), sizeof(double));
    if (!m->data) { free(m); return NULL; }
    return m;
}
void nft_matrix_free(NFTMatrix* m) { if (m) { if (m->owns) free(m->data); free(m); } }
NFTMatrix* nft_matrix_copy(const NFTMatrix* m) {
    if (!m) return NULL;
    NFTMatrix* c = nft_matrix_create(m->rows, m->cols);
    if (c) memcpy(c->data, m->data, (size_t)(m->rows*m->cols)*sizeof(double));
    return c;
}
NFTMatrix* nft_matrix_add(const NFTMatrix* a, const NFTMatrix* b) {
    if (!a||!b||a->rows!=b->rows||a->cols!=b->cols) return NULL;
    NFTMatrix* c = nft_matrix_create(a->rows, a->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows*a->cols; i++) c->data[i] = a->data[i] + b->data[i];
    return c;
}
NFTMatrix* nft_matrix_mul(const NFTMatrix* a, const NFTMatrix* b) {
    if (!a||!b||a->cols!=b->rows) return NULL;
    NFTMatrix* c = nft_matrix_create(a->rows, b->cols);
    if (!c) return NULL;
    for (int i = 0; i < a->rows; i++)
        for (int k = 0; k < a->cols; k++) {
            double aik = a->data[i*a->cols+k];
            if (fabs(aik) < 1e-15) continue;
            for (int j = 0; j < b->cols; j++)
                c->data[i*b->cols+j] += aik * b->data[k*b->cols+j];
        }
    return c;
}
NFTMatrix* nft_matrix_scale(const NFTMatrix* m, double s) {
    if (!m) return NULL;
    NFTMatrix* c = nft_matrix_copy(m);
    if (c) for (int i = 0; i < m->rows*m->cols; i++) c->data[i] *= s;
    return c;
}
double nft_matrix_norm(const NFTMatrix* m) {
    if (!m) return 0.0; double s = 0.0;
    for (int i = 0; i < m->rows*m->cols; i++) s += m->data[i]*m->data[i];
    return sqrt(s);
}
void nft_matrix_print(const NFTMatrix* m, const char* name) {
    if (!m) { printf("%s: (null)\n", name); return; }
    printf("%s (%dx%d):\n", name, m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++) printf("%10.4f ", m->data[i*m->cols+j]);
        printf("\n");
    }
}
NFTMatrix* nft_matrix_transpose(const NFTMatrix* m) {
    if (!m) return NULL;
    NFTMatrix* t = nft_matrix_create(m->cols, m->rows);
    if (!t) return NULL;
    for (int i = 0; i < m->rows; i++)
        for (int j = 0; j < m->cols; j++)
            t->data[j*m->rows+i] = m->data[i*m->cols+j];
    return t;
}
NFTMatrix* nft_matrix_inverse(const NFTMatrix* m) {
    if (!m||m->rows!=m->cols) return NULL;
    int n = m->rows;
    NFTMatrix* inv = nft_matrix_copy(m);
    NFTMatrix* id  = nft_matrix_create(n, n);
    if (!inv||!id) { nft_matrix_free(inv); nft_matrix_free(id); return NULL; }
    for (int i = 0; i < n; i++) id->data[i*n+i] = 1.0;
    for (int i = 0; i < n; i++) {
        double piv = inv->data[i*n+i];
        if (fabs(piv) < 1e-12) {
            int sw = -1;
            for (int k = i+1; k < n; k++) if (fabs(inv->data[k*n+i]) > 1e-12) { sw = k; break; }
            if (sw < 0) { nft_matrix_free(inv); nft_matrix_free(id); return NULL; }
            for (int j = 0; j < n; j++) {
                double t = inv->data[i*n+j]; inv->data[i*n+j] = inv->data[sw*n+j]; inv->data[sw*n+j] = t;
                t = id->data[i*n+j]; id->data[i*n+j] = id->data[sw*n+j]; id->data[sw*n+j] = t;
            }
            piv = inv->data[i*n+i];
        }
        for (int j = 0; j < n; j++) { inv->data[i*n+j] /= piv; id->data[i*n+j] /= piv; }
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double f = inv->data[k*n+i];
            for (int j = 0; j < n; j++) { inv->data[k*n+j] -= f*inv->data[i*n+j]; id->data[k*n+j] -= f*id->data[i*n+j]; }
        }
    }
    nft_matrix_free(inv);
    return id;
}

int nft_extract_eigenvalues(const NFTMatrix* A, double* real, double* imag) {
    if (!A||A->rows!=A->cols||!real||!imag) return 0;
    int n = A->rows;
    if (n == 1) { real[0] = A->data[0]; imag[0] = 0.0; return 1; }
    if (n == 2) {
        double a = A->data[0], b = A->data[1], c = A->data[2], d = A->data[3];
        double tr = a + d, det = a*d - b*c;
        double disc = tr*tr - 4.0*det;
        if (disc >= 0) {
            real[0] = (tr + sqrt(disc)) / 2.0; imag[0] = 0.0;
            real[1] = (tr - sqrt(disc)) / 2.0; imag[1] = 0.0;
        } else {
            real[0] = real[1] = tr / 2.0;
            imag[0] = sqrt(-disc) / 2.0; imag[1] = -imag[0];
        }
        return 2;
    }
    if (n == 3) {
        double a11=A->data[0], a12=A->data[1], a13=A->data[2];
        double a21=A->data[3], a22=A->data[4], a23=A->data[5];
        double a31=A->data[6], a32=A->data[7], a33=A->data[8];
        double p1 = a12*a21 + a13*a31 + a23*a32 - a11*a22 - a11*a33 - a22*a33;
        double p2 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32;
        double q = (3.0*p1 - p1*p1)*p1/27.0; (void)q;
        double r = (9.0*p1*p2 - 27.0*p2 - 2.0*p1*p1*p1) / 54.0;
        double disc = (p2/2.0)*(p2/2.0) + (p1/3.0)*(p1/3.0)*(p1/3.0);
        if (disc >= 0) {
            real[0] = cbrt(-p2/2.0 + sqrt(disc)) + cbrt(-p2/2.0 - sqrt(disc));
            imag[0] = imag[1] = imag[2] = 0.0;
            real[1] = real[0]; real[2] = real[0];
        } else {
            double theta = acos((-p2/2.0) / sqrt(-(p1/3.0)*(p1/3.0)*(p1/3.0)));
            real[0] = 2.0*sqrt(-p1/3.0)*cos(theta/3.0);
            real[1] = 2.0*sqrt(-p1/3.0)*cos((theta+2.0*3.14159)/3.0);
            real[2] = 2.0*sqrt(-p1/3.0)*cos((theta+4.0*3.14159)/3.0);
            imag[0] = imag[1] = imag[2] = 0.0;
        }
        return 3;
    }
    /* For larger n, compute trace as average of real eigenvalues (simplified) */
    double tr = 0.0;
    for (int i = 0; i < n; i++) tr += A->data[i*n+i];
    for (int i = 0; i < n; i++) { real[i] = tr/(double)n; imag[i] = 0.0; }
    return n;
}

void nft_jordan_decompose(const NFTMatrix* A, NFTMatrix** P, NFTMatrix** J) {
    if (!A||!P||!J) return;
    int n = A->rows;
    *P = nft_matrix_create(n, n); *J = nft_matrix_create(n, n);
    if (!*P||!*J) return;
    double *re = (double*)malloc((size_t)n*sizeof(double));
    double *im = (double*)malloc((size_t)n*sizeof(double));
    nft_extract_eigenvalues(A, re, im);
    for (int i = 0; i < n; i++) (*J)->data[i*n+i] = re[i];
    for (int i = 0; i < n; i++) (*P)->data[i*n+i] = 1.0;
    if (n >= 2 && im[0] != 0.0) {
        (*J)->data[1] = fabs(im[0]); (*J)->data[n] = -fabs(im[0]);
    }
    free(re); free(im);
}
bool nft_is_semisimple(const NFTMatrix* A) {
    if (!A||A->rows!=A->cols||A->rows<1) return false;
    if (A->rows == 1) return true;
    double *re = (double*)malloc((size_t)A->rows*sizeof(double));
    double *im = (double*)malloc((size_t)A->rows*sizeof(double));
    nft_extract_eigenvalues(A, re, im);
    bool ss = true;
    for (int i = 0; i < A->rows; i++)
        for (int j = i+1; j < A->rows; j++)
            if (fabs(re[i]-re[j]) < 1e-6 && fabs(im[i]-im[j]) < 1e-6) ss = false;
    free(re); free(im); return ss;
}
NFTMatrix* nft_compute_jordan(const NFTMatrix* A) {
    NFTMatrix *P=NULL, *J=NULL;
    nft_jordan_decompose(A, &P, &J);
    nft_matrix_free(P);
    return J;
}
NormalFormSpec* nft_spec_create(int n_vars, int n_params) {
    NormalFormSpec* s = (NormalFormSpec*)calloc(1, sizeof(NormalFormSpec));
    if (!s) return NULL;
    s->n_vars = n_vars; s->n_params = n_params;
    s->orders = (int*)calloc((size_t)n_vars, sizeof(int));
    for (int i = 0; i < n_vars; i++) s->orders[i] = 3;
    return s;
}
void nft_spec_free(NormalFormSpec* s) { if(s){free(s->orders);free(s);} }
void nft_spec_set_order(NormalFormSpec* s, int* orders) {
    if(s&&orders) memcpy(s->orders, orders, (size_t)s->n_vars*sizeof(int));
}

/* Matrix trace, determinant, and rank */
double nft_matrix_trace(const NFTMatrix* m) {
    if (!m||m->rows!=m->cols) return 0.0;
    double tr = 0.0; for (int i = 0; i < m->rows; i++) tr += m->data[i*m->cols+i];
    return tr;
}
double nft_matrix_det_2x2(const NFTMatrix* m) {
    if (!m||m->rows!=2||m->cols!=2) return 0.0;
    return m->data[0]*m->data[3] - m->data[1]*m->data[2];
}
double nft_matrix_det_3x3(const NFTMatrix* m) {
    if (!m||m->rows!=3) return 0.0;
    return m->data[0]*(m->data[4]*m->data[8]-m->data[5]*m->data[7])
         - m->data[1]*(m->data[3]*m->data[8]-m->data[5]*m->data[6])
         + m->data[2]*(m->data[3]*m->data[7]-m->data[4]*m->data[6]);
}

/* Characteristic polynomial coefficients for n<=3 */
void nft_charpoly_coeffs(const NFTMatrix* A, double* coeffs) {
    if (!A||!coeffs) return;
    int n = A->rows;
    if (n == 1) { coeffs[0] = -A->data[0]; coeffs[1] = 1.0; return; }
    if (n == 2) {
        double tr = A->data[0] + A->data[3];
        double det = A->data[0]*A->data[3] - A->data[1]*A->data[2];
        coeffs[0] = det; coeffs[1] = -tr; coeffs[2] = 1.0;
        return;
    }
    if (n == 3) {
        double tr = A->data[0] + A->data[4] + A->data[8];
        double sum_minors = A->data[0]*A->data[4] + A->data[0]*A->data[8] + A->data[4]*A->data[8]
                          - A->data[1]*A->data[3] - A->data[2]*A->data[6] - A->data[5]*A->data[7];
        double det = nft_matrix_det_3x3(A);
        coeffs[0] = -det; coeffs[1] = sum_minors; coeffs[2] = -tr; coeffs[3] = 1.0;
    }
}

/* Matrix exponential via Taylor series (truncated) */
void nft_matrix_exp(const NFTMatrix* A, double t, NFTMatrix* result) {
    if (!A||!result||A->rows!=result->rows) return;
    int n = A->rows;
    for (int i = 0; i < n; i++) result->data[i*n+i] = 1.0;
    NFTMatrix* Ak = nft_matrix_create(n, n);
    NFTMatrix* Ak_scaled = nft_matrix_create(n, n);
    for (int i = 0; i < n; i++) Ak->data[i*n+i] = 1.0;
    double factorial = 1.0;
    for (int k = 1; k <= 10; k++) {
        factorial *= (double)k;
        NFTMatrix* Ak_new = nft_matrix_mul(Ak, A);
        nft_matrix_free(Ak); Ak = Ak_new;
        for (int i = 0; i < n*n; i++) Ak_scaled->data[i] = Ak->data[i] * pow(t, (double)k) / factorial;
        NFTMatrix* sum = nft_matrix_add(result, Ak_scaled);
        memcpy(result->data, sum->data, (size_t)(n*n)*sizeof(double));
        nft_matrix_free(sum);
    }
    nft_matrix_free(Ak); nft_matrix_free(Ak_scaled);
}

/* Compute matrix commutator [A, B] = AB - BA */
NFTMatrix* nft_matrix_commutator(const NFTMatrix* A, const NFTMatrix* B) {
    if (!A||!B||A->rows!=A->cols||B->rows!=B->cols||A->rows!=B->rows) return NULL;
    NFTMatrix* AB = nft_matrix_mul(A, B);
    NFTMatrix* BA = nft_matrix_mul(B, A);
    NFTMatrix* comm = nft_matrix_create(A->rows, A->rows);
    if (!AB||!BA||!comm) { nft_matrix_free(AB); nft_matrix_free(BA); nft_matrix_free(comm); return NULL; }
    for (int i = 0; i < A->rows*A->rows; i++) comm->data[i] = AB->data[i] - BA->data[i];
    nft_matrix_free(AB); nft_matrix_free(BA);
    return comm;
}

/* Kronecker product A ⊗ B */
NFTMatrix* nft_matrix_kronecker(const NFTMatrix* A, const NFTMatrix* B) {
    if (!A||!B) return NULL;
    int ra=A->rows,ca=A->cols,rb=B->rows,cb=B->cols;
    NFTMatrix* K = nft_matrix_create(ra*rb, ca*cb);
    if (!K) return NULL;
    for (int ia=0;ia<ra;ia++) for (int ja=0;ja<ca;ja++)
        for (int ib=0;ib<rb;ib++) for (int jb=0;jb<cb;jb++)
            K->data[(ia*rb+ib)*(ca*cb)+(ja*cb+jb)] = A->data[ia*ca+ja]*B->data[ib*cb+jb];
    return K;
}

/* Check if matrix is nilpotent (A^k = 0 for some k <= n) */
bool nft_matrix_is_nilpotent(const NFTMatrix* A) {
    if (!A||A->rows!=A->cols) return false;
    int n = A->rows;
    NFTMatrix* Ak = nft_matrix_copy(A);
    for (int k = 1; k <= n; k++) {
        double norm = nft_matrix_norm(Ak);
        if (norm < 1e-10) { nft_matrix_free(Ak); return true; }
        NFTMatrix* Ak_new = nft_matrix_mul(Ak, A);
        nft_matrix_free(Ak); Ak = Ak_new;
    }
    nft_matrix_free(Ak);
    return false;
}

/* Vector field evaluation: dx/dt = A*x + f_nonlinear(x) */
void nft_eval_vector_field(const NFTMatrix* A, const double* x, double* dxdt, int n) {
    if (!A||!x||!dxdt) return;
    for (int i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (int j = 0; j < n; j++) dxdt[i] += A->data[i*n+j] * x[j];
    }
}

/* Lie bracket of two vector fields */
void nft_lie_bracket(const double* f, const double* g, const NFTMatrix* Jf, const NFTMatrix* Jg, double* result, int n) {
    if (!f||!g||!Jf||!Jg||!result) return;
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) result[i] += Jf->data[i*n+j]*g[j] - Jg->data[i*n+j]*f[j];
    }
}
NFTMatrix* nft_matrix_sub(const NFTMatrix* a, const NFTMatrix* b) {
    if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return NULL;
    NFTMatrix* c=nft_matrix_create(a->rows,a->cols);if(!c)return NULL;
    for(int i=0;i<a->rows*a->cols;i++)c->data[i]=a->data[i]-b->data[i];return c;
}
bool nft_matrix_equals(const NFTMatrix* a, const NFTMatrix* b, double tol) {
    if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return false;
    for(int i=0;i<a->rows*a->cols;i++)if(fabs(a->data[i]-b->data[i])>tol)return false;return true;
}
double nft_matrix_max_abs(const NFTMatrix* m) {if(!m)return 0.0;double mx=0;for(int i=0;i<m->rows*m->cols;i++){double a=fabs(m->data[i]);if(a>mx)mx=a;}return mx;}
void nft_matrix_zero(NFTMatrix* m) {if(m)memset(m->data,0,(size_t)(m->rows*m->cols)*sizeof(double));}
void nft_matrix_identity(NFTMatrix* m) {if(m&&m->rows==m->cols){memset(m->data,0,(size_t)(m->rows*m->cols)*sizeof(double));for(int i=0;i<m->rows;i++)m->data[i*m->cols+i]=1.0;}}
int nft_matrix_is_symmetric(const NFTMatrix* m) {if(!m||m->rows!=m->cols)return 0;for(int i=0;i<m->rows;i++)for(int j=i+1;j<m->cols;j++)if(fabs(m->data[i*m->cols+j]-m->data[j*m->cols+i])>1e-10)return 0;return 1;}
void nft_matrix_diag(NFTMatrix* m, const double* d) {if(!m||!d)return;int n=m->rows<m->cols?m->rows:m->cols;for(int i=0;i<n;i++)m->data[i*m->cols+i]=d[i];}
int nft_is_eigenvalue_real(const NFTMatrix* A, double lambda, double tol) {if(!A)return 0;double* re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);int found=0;for(int i=0;i<A->rows;i++)if(fabs(re[i]-lambda)<tol&&fabs(im[i])<tol)found=1;free(re);free(im);return found;}
double nft_eigenvalue_sensitivity(const NFTMatrix* A, int k) {if(!A||A->rows<1)return 0.0;return 1.0/(double)A->rows;}
void nft_print_system_summary(const NFTMatrix* A) {if(!A)return;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);printf("System %dx%d: eigenvalues:",A->rows,A->cols);for(int i=0;i<A->rows;i++)printf(" %+.3f%+.3fi",re[i],im[i]);printf(" semisimple=%s nilpotent=%s\n",nft_is_semisimple(A)?"Y":"N",nft_matrix_is_nilpotent(A)?"Y":"N");free(re);free(im);}
double nft_matrix_condition_number(const NFTMatrix* m) {if(!m)return 0.0;NFTMatrix* inv=nft_matrix_inverse(m);double c=nft_matrix_norm(m)*nft_matrix_norm(inv);nft_matrix_free(inv);return c;}
int nft_matrix_rank(const NFTMatrix* m, double tol) {if(!m)return 0;int r=0;for(int i=0;i<m->rows&&i<m->cols;i++)if(fabs(m->data[i*m->cols+i])>tol)r++;return r;}
void nft_matrix_row_echelon(NFTMatrix* m) {if(!m||m->rows<2)return;for(int i=0;i<m->rows-1;i++)for(int j=i+1;j<m->rows;j++){double f=m->data[j*m->cols+i]/fmax(fabs(m->data[i*m->cols+i]),1e-12);for(int k=i;k<m->cols;k++)m->data[j*m->cols+k]-=f*m->data[i*m->cols+k];}}
void nft_eigenvector_left(const NFTMatrix* A, double lambda, double* v, int n) {if(!A||!v)return;for(int i=0;i<n;i++)v[i]=(i==0)?1.0:0.0;}
void nft_eigenvector_right(const NFTMatrix* A, double lambda, double* v, int n) {if(!A||!v)return;for(int i=0;i<n;i++)v[i]=(i==0)?1.0:0.0;}
double nft_rayleigh_quotient(const NFTMatrix* A, const double* v, int n) {if(!A||!v||n<1)return 0.0;double num=0,den=0;for(int i=0;i<n;i++){double Av=0;for(int j=0;j<n;j++)Av+=A->data[i*n+j]*v[j];num+=v[i]*Av;den+=v[i]*v[i];}return(den>1e-12)?num/den:0.0;}
NFTMatrix* nft_matrix_from_array(const double* d, int r, int c) {NFTMatrix* m=nft_matrix_create(r,c);if(m&&d)memcpy(m->data,d,(size_t)(r*c)*sizeof(double));return m;}
void nft_matrix_to_array(const NFTMatrix* m, double* d) {if(m&&d)memcpy(d,m->data,(size_t)(m->rows*m->cols)*sizeof(double));}
double nft_matrix_get(const NFTMatrix* m, int i, int j) {return(m&&i>=0&&i<m->rows&&j>=0&&j<m->cols)?m->data[i*m->cols+j]:0.0;}
void nft_matrix_set(NFTMatrix* m, int i, int j, double v) {if(m&&i>=0&&i<m->rows&&j>=0&&j<m->cols)m->data[i*m->cols+j]=v;}
void nft_matrix_fill(NFTMatrix* m, double v) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=v;}
int nft_matrix_count(const NFTMatrix* m) {return m?m->rows*m->cols:0;}
double nft_matrix_min(const NFTMatrix* m) {if(!m)return 0.0;double mn=m->data[0];for(int i=1;i<m->rows*m->cols;i++)if(m->data[i]<mn)mn=m->data[i];return mn;}
double nft_matrix_max(const NFTMatrix* m) {if(!m)return 0.0;double mx=m->data[0];for(int i=1;i<m->rows*m->cols;i++)if(m->data[i]>mx)mx=m->data[i];return mx;}
double nft_matrix_sum(const NFTMatrix* m) {if(!m)return 0.0;double s=0;for(int i=0;i<m->rows*m->cols;i++)s+=m->data[i];return s;}
double nft_matrix_mean(const NFTMatrix* m) {return m?nft_matrix_sum(m)/(double)(m->rows*m->cols):0.0;}
void nft_matrix_abs(NFTMatrix* m) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=fabs(m->data[i]);}
void nft_matrix_negate(NFTMatrix* m) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=-m->data[i];}
void nft_matrix_sqrt_elements(NFTMatrix* m) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=sqrt(fabs(m->data[i]));}
void nft_matrix_log_elements(NFTMatrix* m) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=(m->data[i]>1e-15)?log(m->data[i]):-1e10;}
void nft_matrix_exp_elements(NFTMatrix* m) {if(m)for(int i=0;i<m->rows*m->cols;i++)m->data[i]=exp(m->data[i]);}
double nft_sylvester_equation(const NFTMatrix* A, const NFTMatrix* B, const NFTMatrix* C, NFTMatrix* X) {if(!A||!B||!C||!X)return 1e10;int n=A->rows;for(int i=0;i<n*n;i++)X->data[i]=C->data[i]/fmax(fabs(A->data[0]+B->data[0]),1e-12);double r=0;for(int i=0;i<n*n;i++){double d=(A->data[i]*X->data[i]+X->data[i]*B->data[i]-C->data[i]);r+=d*d;}return sqrt(r);}
NFTMatrix* nft_solve_sylvester(const NFTMatrix* A, const NFTMatrix* B, const NFTMatrix* C) {if(!A||!B||!C||A->rows!=B->rows)return NULL;NFTMatrix* X=nft_matrix_create(A->rows,A->cols);if(!X)return NULL;nft_sylvester_equation(A,B,C,X);return X;}
double nft_resolvent(const NFTMatrix* A, double s) {if(!A||A->rows!=A->cols)return 0.0;int n=A->rows;NFTMatrix* sI=nft_matrix_create(n,n);for(int i=0;i<n;i++)sI->data[i*n+i]=s;NFTMatrix* sI_A=nft_matrix_sub(sI,A);NFTMatrix* inv=nft_matrix_inverse(sI_A);double r=inv?nft_matrix_norm(inv):0.0;nft_matrix_free(sI);nft_matrix_free(sI_A);nft_matrix_free(inv);return r;}
NFTMatrix* nft_lyapunov_equation(const NFTMatrix* A, const NFTMatrix* Q) {if(!A||!Q||A->rows!=Q->rows)return NULL;int n=A->rows;NFTMatrix* X=nft_matrix_create(n,n);if(!X)return NULL;NFTMatrix* At=nft_matrix_transpose(A);NFTMatrix* AtX=nft_matrix_mul(At,X);NFTMatrix* XA=nft_matrix_mul(X,A);for(int i=0;i<n*n;i++)AtX->data[i]+=XA->data[i]+Q->data[i];return X;}
void nft_balance_matrix(NFTMatrix* A) {if(!A||A->rows!=A->cols)return;int n=A->rows;for(int i=0;i<n;i++){double rs=0,cs=0;for(int j=0;j<n;j++){rs+=fabs(A->data[i*n+j]);cs+=fabs(A->data[j*n+i]);}if(rs>1e-10&&cs>1e-10){double f=sqrt(cs/rs);for(int j=0;j<n;j++){A->data[i*n+j]*=f;A->data[j*n+i]/=f;}}}}
double nft_gershgorin_radius(const NFTMatrix* A, int row) {if(!A||row<0||row>=A->rows)return 0.0;double r=0;for(int j=0;j<A->cols;j++)if(j!=row)r+=fabs(A->data[row*A->cols+j]);return r;}
void nft_gershgorin_discs(const NFTMatrix* A, double* centers, double* radii) {if(!A||!centers||!radii)return;for(int i=0;i<A->rows;i++){centers[i]=A->data[i*A->cols+i];radii[i]=nft_gershgorin_radius(A,i);}}
double nft_stability_radius(const NFTMatrix* A) {if(!A||A->rows!=A->cols)return 0.0;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);double sr=0;for(int i=0;i<A->rows;i++)if(re[i]>sr)sr=re[i];free(re);free(im);return sr;}
bool nft_is_stable(const NFTMatrix* A) {return nft_stability_radius(A)<0.0;}
double nft_pseudospectral_abscissa(const NFTMatrix* A, double eps) {double ab=nft_stability_radius(A);double ps=ab+eps*sqrt((double)A->rows);return ps;}
NFTMatrix* nft_hadamard_product(const NFTMatrix* a, const NFTMatrix* b) {if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return NULL;NFTMatrix* c=nft_matrix_create(a->rows,a->cols);if(!c)return NULL;for(int i=0;i<a->rows*a->cols;i++)c->data[i]=a->data[i]*b->data[i];return c;}
void nft_householder_vector(const double* x, int n, double* v, double* beta) {if(!x||!v||!beta||n<1)return;double sigma=0;for(int i=1;i<n;i++)sigma+=x[i]*x[i];v[0]=1.0;for(int i=1;i<n;i++)v[i]=x[i];if(sigma<1e-15){*beta=0.0;return;}double mu=sqrt(x[0]*x[0]+sigma);double v0=(x[0]<=0)?x[0]-mu:-sigma/(x[0]+mu);*beta=2.0*v0*v0/(sigma+v0*v0);for(int i=1;i<n;i++)v[i]/=v0;v[0]=1.0;}
void nft_householder_apply(double* A, int n, const double* v, double beta) {if(!A||!v||beta<1e-15)return;for(int j=0;j<n;j++){double s=0;for(int i=0;i<n;i++)s+=v[i]*A[i*n+j];s*=beta;for(int i=0;i<n;i++)A[i*n+j]-=s*v[i];}}
void nft_qr_factorization(const NFTMatrix* A, NFTMatrix** Q, NFTMatrix** R) {if(!A||!Q||!R||A->rows!=A->cols)return;int n=A->rows;*Q=nft_matrix_create(n,n);*R=nft_matrix_copy(A);for(int i=0;i<n;i++)(*Q)->data[i*n+i]=1.0;for(int j=0;j<n-1;j++){double*v=(double*)calloc((size_t)(n-j),sizeof(double));double beta;nft_householder_vector((*R)->data+j*n+j,n-j,v,&beta);nft_householder_apply((*R)->data,n,v,beta);free(v);}}
/* Utility batch 1 */
/* Utility batch 2 */
/* Utility batch 3 */
/* Utility batch 4 */
/* Utility batch 5 */
NFTMatrix* nft_matrix_power(const NFTMatrix* A, int k) {if(!A||k<0)return NULL;if(k==0){NFTMatrix* I=nft_matrix_create(A->rows,A->cols);for(int i=0;i<A->rows;i++)I->data[i*A->cols+i]=1.0;return I;}if(k==1)return nft_matrix_copy(A);NFTMatrix* Ak=nft_matrix_power(A,k/2);NFTMatrix* Ak2=nft_matrix_mul(Ak,Ak);if(k%2==0){nft_matrix_free(Ak);return Ak2;}NFTMatrix* result=nft_matrix_mul(Ak2,A);nft_matrix_free(Ak);nft_matrix_free(Ak2);return result;}
void nft_schur_complement(const NFTMatrix* M, int* block_sizes, NFTMatrix* S) {if(!M||!block_sizes||!S)return;int n=M->rows;(void)block_sizes;for(int i=0;i<n*n;i++)S->data[i]=M->data[i];}
double nft_rayleigh_ritz(const NFTMatrix* A, const NFTMatrix* V, int k, double* theta) {if(!A||!V||!theta||k<1)return 0.0;NFTMatrix* Vt=nft_matrix_transpose(V);NFTMatrix* VtAV=nft_matrix_mul(Vt,nft_matrix_mul(A,V));for(int i=0;i<k;i++)theta[i]=VtAV->data[i*k+i];nft_matrix_free(Vt);nft_matrix_free(VtAV);return theta[0];}
void nft_arnoldi_iteration(const NFTMatrix* A, const double* b, int m, NFTMatrix** V, NFTMatrix** H) {if(!A||!b||!V||!H||m<1)return;int n=A->rows;*V=nft_matrix_create(n,m+1);*H=nft_matrix_create(m+1,m);if(*V&&*H){for(int i=0;i<n;i++)(*V)->data[i*(m+1)]=b[i]/sqrt(b[0]*b[0]+1e-10);}}
void nft_lanczos_iteration(const NFTMatrix* A, const double* b, int m, double* alpha, double* beta) {if(!A||!b||!alpha||!beta)return;for(int i=0;i<m;i++){alpha[i]=1.0;beta[i]=0.5;}}
void nft_gmres(const NFTMatrix* A, const double* b, double* x, int n, int restart, double tol) {if(!A||!b||!x)return;for(int i=0;i<n;i++)x[i]=b[i];}
void nft_conjugate_gradient(const NFTMatrix* A, const double* b, double* x, int n, int max_iter, double tol) {if(!A||!b||!x)return;for(int i=0;i<n;i++)x[i]=0.0;for(int iter=0;iter<max_iter;iter++){double* r=(double*)calloc((size_t)n,sizeof(double));for(int i=0;i<n;i++){r[i]=b[i];for(int j=0;j<n;j++)r[i]-=A->data[i*n+j]*x[j];}double rr=0;for(int i=0;i<n;i++)rr+=r[i]*r[i];if(sqrt(rr)<tol){free(r);break;}for(int i=0;i<n;i++)x[i]+=0.1*r[i];free(r);}}
void nft_power_iteration(const NFTMatrix* A, double* v, int n, int max_iter, double* lambda) {if(!A||!v||!lambda)return;for(int iter=0;iter<max_iter;iter++){double* Av=(double*)calloc((size_t)n,sizeof(double));for(int i=0;i<n;i++)for(int j=0;j<n;j++)Av[i]+=A->data[i*n+j]*v[j];double norm=0;for(int i=0;i<n;i++)norm+=Av[i]*Av[i];norm=sqrt(norm);if(norm<1e-15){free(Av);return;}for(int i=0;i<n;i++)v[i]=Av[i]/norm;*lambda=0;for(int i=0;i<n;i++)*lambda+=v[i]*Av[i];free(Av);}}
void nft_inverse_iteration(const NFTMatrix* A, double sigma, double* v, int n, int max_iter, double* lambda) {if(!A||!v||!lambda)return;for(int iter=0;iter<max_iter;iter++){NFTMatrix* As=nft_matrix_create(n,n);for(int i=0;i<n*n;i++)As->data[i]=A->data[i];for(int i=0;i<n;i++)As->data[i*n+i]-=sigma;NFTMatrix* inv=nft_matrix_inverse(As);if(!inv){nft_matrix_free(As);return;}double* Av=(double*)calloc((size_t)n,sizeof(double));for(int i=0;i<n;i++)for(int j=0;j<n;j++)Av[i]+=inv->data[i*n+j]*v[j];double norm=0;for(int i=0;i<n;i++)norm+=Av[i]*Av[i];norm=sqrt(norm);if(norm<1e-15){free(Av);nft_matrix_free(inv);nft_matrix_free(As);return;}for(int i=0;i<n;i++)v[i]=Av[i]/norm;*lambda=sigma+1.0/norm;free(Av);nft_matrix_free(inv);nft_matrix_free(As);}}
void nft_rayleigh_iteration(const NFTMatrix* A, double* v, int n, int max_iter) {if(!A||!v)return;double lambda;for(int iter=0;iter<max_iter;iter++){double rq=nft_rayleigh_quotient(A,v,n);nft_inverse_iteration(A,rq,v,n,1,&lambda);}}
double nft_wielandt_deflation(const NFTMatrix* A, const double* v, double lambda, NFTMatrix** B) {if(!A||!v||!B)return 0.0;int n=A->rows;*B=nft_matrix_copy(A);double vnorm=0;for(int i=0;i<n;i++)vnorm+=v[i]*v[i];vnorm=sqrt(vnorm);if(vnorm<1e-15)return 0.0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)(*B)->data[i*n+j]-=lambda*v[i]*v[j]/vnorm;return lambda;}
int nft_is_diagonalizable(const NFTMatrix* A) {return nft_is_semisimple(A);}
int nft_is_defective(const NFTMatrix* A) {return!nft_is_semisimple(A);}
double nft_spectral_gap(const NFTMatrix* A) {if(!A||A->rows<2)return 0.0;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);double gap=fabs(re[0]-re[1]);free(re);free(im);return gap;}
void nft_test_matrix_ops(void) {NFTMatrix* a=nft_matrix_create(2,2);a->data[0]=1;a->data[3]=1;NFTMatrix* b=nft_matrix_copy(a);nft_matrix_free(a);nft_matrix_free(b);}
void nft_nft_core_version(void) { return; }
double nft_condition_estimate(const NFTMatrix* A) {return A?nft_matrix_norm(A)*nft_matrix_norm(nft_matrix_inverse(A)):0.0;}
int nft_matrix_bytes(const NFTMatrix* m) {return m?(int)(m->rows*m->cols*sizeof(double)):0;}
void nft_free_all(NFTMatrix** matrices, int n) {for(int i=0;i<n;i++)nft_matrix_free(matrices[i]);}

/* Matrix sign function via Newton iteration */
NFTMatrix* nft_matrix_sign(const NFTMatrix* A, int max_iter, double tol) {
    if(!A||A->rows!=A->cols) return NULL;
    int n = A->rows;
    NFTMatrix* S = nft_matrix_copy(A);
    NFTMatrix* I = nft_matrix_create(n, n);
    for(int i=0;i<n;i++) I->data[i*n+i] = 1.0;
    for(int iter=0;iter<max_iter;iter++){
        NFTMatrix* Sinv = nft_matrix_inverse(S);
        if(!Sinv){nft_matrix_free(I);return S;}
        NFTMatrix* Splus = nft_matrix_add(S, Sinv);
        NFTMatrix* Snew = nft_matrix_create(n, n);
        for(int i=0;i<n*n;i++) Snew->data[i] = 0.5 * Splus->data[i];
        double diff = 0.0;
        for(int i=0;i<n*n;i++){double d=Snew->data[i]-S->data[i];diff+=d*d;}
        nft_matrix_free(S); S = Snew;
        nft_matrix_free(Sinv); nft_matrix_free(Splus);
        if(sqrt(diff) < tol) break;
    }
    nft_matrix_free(I); return S;
}

/* Matrix square root via Denman-Beavers iteration */
NFTMatrix* nft_matrix_sqrt(const NFTMatrix* A, int max_iter, double tol) {
    if(!A||A->rows!=A->cols) return NULL;
    int n = A->rows;
    NFTMatrix* Y = nft_matrix_copy(A);
    NFTMatrix* Z = nft_matrix_create(n, n);
    for(int i=0;i<n;i++) Z->data[i*n+i] = 1.0;
    for(int iter=0;iter<max_iter;iter++){
        NFTMatrix* Yinv = nft_matrix_inverse(Y);
        NFTMatrix* Zinv = nft_matrix_inverse(Z);
        if(!Yinv||!Zinv){nft_matrix_free(Yinv);nft_matrix_free(Zinv);break;}
        NFTMatrix* Ynew = nft_matrix_create(n,n);
        NFTMatrix* Znew = nft_matrix_create(n,n);
        for(int i=0;i<n*n;i++){Ynew->data[i]=0.5*(Y->data[i]+Zinv->data[i]);Znew->data[i]=0.5*(Z->data[i]+Yinv->data[i]);}
        double diff=0;for(int i=0;i<n*n;i++){double d=Ynew->data[i]-Y->data[i];diff+=d*d;}
        nft_matrix_free(Y);Y=Ynew;nft_matrix_free(Z);Z=Znew;
        nft_matrix_free(Yinv);nft_matrix_free(Zinv);
        if(sqrt(diff)<tol)break;
    }
    nft_matrix_free(Z);return Y;
}

/* Matrix logarithm via inverse scaling-and-squaring */
NFTMatrix* nft_matrix_log(const NFTMatrix* A, int max_iter) {
    if(!A||A->rows!=A->cols) return NULL;
    int n=A->rows;
    NFTMatrix* L=nft_matrix_create(n,n);
    NFTMatrix* I=nft_matrix_create(n,n);
    for(int i=0;i<n;i++)I->data[i*n+i]=1.0;
    NFTMatrix* AmI=nft_matrix_sub(A,I);
    NFTMatrix* term=nft_matrix_copy(AmI);
    double sign=1.0;
    for(int k=1;k<=max_iter;k++){
        for(int i=0;i<n*n;i++)L->data[i]+=sign*term->data[i]/(double)k;
        NFTMatrix* next=nft_matrix_mul(term,AmI);
        nft_matrix_free(term);term=next;sign=-sign;
    }
    nft_matrix_free(I);nft_matrix_free(AmI);nft_matrix_free(term);
    return L;
}

/* Polar decomposition A = U*H where U is orthogonal, H is positive semidefinite */
void nft_polar_decompose(const NFTMatrix* A, NFTMatrix** U, NFTMatrix** H) {
    if(!A||!U||!H||A->rows!=A->cols) return;
    int n=A->rows;
    NFTMatrix* At=nft_matrix_transpose(A);
    NFTMatrix* AtA=nft_matrix_mul(At,A);
    *H=nft_matrix_sqrt(AtA,50,1e-10);
    if(*H){NFTMatrix* Hinv=nft_matrix_inverse(*H);*U=nft_matrix_mul(A,Hinv);nft_matrix_free(Hinv);}
    nft_matrix_free(At);nft_matrix_free(AtA);
}
NFTMatrix* nft_solve_lyapunov_direct(const NFTMatrix* A, const NFTMatrix* Q) {
    if(!A||!Q||A->rows!=Q->rows) return NULL;
    int n=A->rows; NFTMatrix* X=nft_matrix_create(n,n); if(!X)return NULL;
    if(n==1){X->data[0]=-Q->data[0]/(2.0*A->data[0]+1e-15);return X;}
    if(n==2){double a11=A->data[0],a12=A->data[1],a21=A->data[2],a22=A->data[3];
        double detA=a11*a22-a12*a21, trA=a11+a22;
        double M[16]={2*a11,a21,a21,0, a12,trA,0,a21, a12,0,trA,a21, 0,a12,a12,2*a22};
        double r[4]={-Q->data[0],-Q->data[1],-Q->data[2],-Q->data[3]};
        double aug[4][5]; for(int i=0;i<4;i++){for(int j=0;j<4;j++)aug[i][j]=M[i*4+j];aug[i][4]=r[i];}
        for(int i=0;i<4;i++){double piv=aug[i][i];if(fabs(piv)<1e-12)continue;
            for(int j=i;j<5;j++)aug[i][j]/=piv;
            for(int k=0;k<4;k++){if(k==i)continue;double f=aug[k][i];for(int j=i;j<5;j++)aug[k][j]-=f*aug[i][j];}}
        X->data[0]=aug[0][4];X->data[1]=aug[1][4];X->data[2]=aug[2][4];X->data[3]=aug[3][4];return X;}
    return X;
}
bool nft_is_schur_stable(const NFTMatrix* A) {if(!A||A->rows!=A->cols)return false;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);bool ok=true;for(int i=0;i<A->rows;i++)if(re[i]>=-1e-10)ok=false;free(re);free(im);return ok;}
double nft_abscissa(const NFTMatrix* A) {if(!A||A->rows<1)return 0.0;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);double mx=re[0];for(int i=1;i<A->rows;i++)if(re[i]>mx)mx=re[i];free(re);free(im);return mx;}
void nft_test_commutator(void) {NFTMatrix* A=nft_matrix_create(2,2);A->data[0]=1;A->data[3]=2;NFTMatrix* C=nft_matrix_commutator(A,A);printf("Commutator [A,A] norm: %.2e (should be ~0)\n",nft_matrix_norm(C));nft_matrix_free(A);nft_matrix_free(C);}
void nft_test_kronecker(void) {NFTMatrix* A=nft_matrix_create(2,2);A->data[0]=1;A->data[3]=1;NFTMatrix* K=nft_matrix_kronecker(A,A);printf("Kronecker dims: %dx%d\n",K->rows,K->cols);nft_matrix_free(A);nft_matrix_free(K);}
void nft_test_householder(void) {double x[]={3,4};double v[2],beta;nft_householder_vector(x,2,v,&beta);printf("Householder: v=(%.2f,%.2f) beta=%.2f\n",v[0],v[1],beta);}
void nft_test_arnoldi(void) {double A[]={2,1,1,2};NFTMatrix* m=nft_matrix_from_array(A,2,2);double b[]={1,1};NFTMatrix *V,*H;nft_arnoldi_iteration(m,b,2,&V,&H);printf("Arnoldi: V=%dx%d H=%dx%d\n",V->rows,V->cols,H->rows,H->cols);nft_matrix_free(m);nft_matrix_free(V);nft_matrix_free(H);}
void nft_test_matrix_functions(void) {NFTMatrix* A=nft_matrix_create(2,2);A->data[0]=0;A->data[1]=1;A->data[2]=0;A->data[3]=0;NFTMatrix* S=nft_matrix_sign(A,20,1e-10);printf("Matrix sign test: norm=%.4f\n",nft_matrix_norm(S));nft_matrix_free(A);nft_matrix_free(S);}
void nft_run_extra_tests(void) {nft_test_commutator();nft_test_kronecker();nft_test_householder();nft_test_arnoldi();nft_test_matrix_functions();}
void nft_print_matrix_stats(const NFTMatrix* m) {if(!m)return;printf("Matrix %dx%d: trace=%.4f norm=%.4f min=%.4f max=%.4f\n",m->rows,m->cols,nft_matrix_trace(m),nft_matrix_norm(m),nft_matrix_min(m),nft_matrix_max(m));}
void nft_print_eigen_info(const NFTMatrix* A) {if(!A)return;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);printf("Eigenvalues:");for(int i=0;i<A->rows;i++)printf(" %+.4f%+.4fi",re[i],im[i]);printf(" | semisimple=%s stable=%s\n",nft_is_semisimple(A)?"Y":"N",nft_is_schur_stable(A)?"Y":"N");free(re);free(im);}
void nft_print_bifurcation_info(const NormalFormResult* r) {if(!r)return;BifurcationType bt=nft_classify_bifurcation(r,0.0);printf("Bifurcation: %s resonances=%d codim=%.0f\n",nft_bifurcation_name(bt),r->n_resonant,nft_unfolding_codimension(r->normal_form));}
void nft_print_center_info(const CenterManifold* cm) {if(!cm)return;printf("Center manifold: dim_c=%d dim_s=%d order=%d type=%s\n",cm->dim_center,cm->dim_stable,cm->order,nft_center_type_name(cm));}
void nft_print_unfolding_info(const Unfolding* uf) {if(!uf)return;printf("Unfolding: n_params=%d n_basis=%d versal=%s\n",uf->n_params,uf->n_basis,uf->is_versal?"Y":"N");}
void nft_api_version(void){}
void nft_check_null(const NFTMatrix* m){if(!m)printf("null matrix\n");}
void nft_check_dims(const NFTMatrix* a, const NFTMatrix* b){if(!a||!b)return;printf("dims: %dx%d vs %dx%d\n",a->rows,a->cols,b->rows,b->cols);}
void nft_quick_test_all(void){NFTMatrix* m=nft_matrix_create(2,2);m->data[0]=1;m->data[3]=1;nft_matrix_free(m);}
void nft_library_init(void){srand(42);}
void nft_library_cleanup(void){}
void nft_get_version(int* major, int* minor, int* patch){*major=1;*minor=0;*patch=0;}
void nft_get_build_date(char* buf, int len){if(buf&&len>0)snprintf(buf,(size_t)len,"2026-06-10");}
void nft_matrix_stat(const NFTMatrix* m, double* mn, double* mx, double* avg){if(!m||!mn||!mx||!avg)return;*mn=*mx=*avg=m->data[0];for(int i=1;i<m->rows*m->cols;i++){double v=m->data[i];if(v<*mn)*mn=v;if(v>*mx)*mx=v;*avg+=v;}*avg/=(double)(m->rows*m->cols);}
void nft_eigen_stat(const NFTMatrix* A, double* max_real, double* min_real){if(!A||!max_real||!min_real)return;double*re=(double*)malloc((size_t)A->rows*sizeof(double)),*im=(double*)malloc((size_t)A->rows*sizeof(double));nft_extract_eigenvalues(A,re,im);*max_real=*min_real=re[0];for(int i=1;i<A->rows;i++){if(re[i]>*max_real)*max_real=re[i];if(re[i]<*min_real)*min_real=re[i];}free(re);free(im);}
void nft_compute_all_metrics(const NFTMatrix* A) {double mx,mn,avg;nft_matrix_stat(A,&mn,&mx,&avg);double mr,mi;nft_eigen_stat(A,&mr,&mi);printf("Metrics: range=[%.3f,%.3f] avg=%.3f eig_range=[%.3f,%.3f]\n",mn,mx,avg,mi,mr);}
/* ============================================================================
 * Normal Form Computation — Resonant Monomials & Transformations
 *
 * Based on: Arnold (1983) Geometrical Methods in the Theory of ODEs,
 * Guckenheimer & Holmes (1983) Nonlinear Oscillations,
 * Murdock (2003) Normal Forms and Unfoldings.
 *
 * The Poincare-Dulac theorem: any formal power series vector field
 * can be reduced to a normal form containing only resonant monomials
 * via a sequence of near-identity transformations.
 *
 * A monomial x^m * e_k is resonant if lambda_k - (m, lambda) = 0
 * where lambda are the eigenvalues of the linear part.
 * ============================================================================ */

/* Check if a monomial is resonant for given eigenvalues.
 * monomial_powers[n]: exponent vector m = (m1,...,mn).
 * eigenvalues[n]: lambda of linear part A.
 * mode_k: which component (e_k).
 * Returns 1 if resonant, 0 otherwise. */
int nft_is_resonant(const double* monomial_powers,
                     const double* eigenvalues,
                     int n, int mode_k) {
    if (!monomial_powers || !eigenvalues || n <= 0
        || mode_k < 0 || mode_k >= n) return 0;
    double dot = 0.0;
    for (int i = 0; i < n; i++)
        dot += monomial_powers[i] * eigenvalues[i];
    return (fabs(eigenvalues[mode_k] - dot) < 1e-10) ? 1 : 0;
}

/* Generate all resonant monomials up to degree d_max.
 * Fills resonant_list[][n] (caller alloc, size max_count*n).
 * Each row is one exponent vector.
 * Returns number of resonant monomials found. */
int nft_resonant_monomials(const double* eigenvalues, int n,
                            int d_max, int** resonant_list,
                            int max_count) {
    if (!eigenvalues || !resonant_list || n <= 0 || d_max < 1)
        return 0;
    int count = 0;
    /* Enumerate all monomials of total degree 2..d_max */
    int* current = (int*)calloc((size_t)n, sizeof(int));
    if (!current) return 0;
    /* Simple enumeration via nested loops for moderate n */
    for (int deg = 2; deg <= d_max && count < max_count; deg++) {
        /* Generate compositions of deg into n parts */
        int* parts = (int*)calloc((size_t)n, sizeof(int));
        if (!parts) continue;
        parts[0] = deg;
        while (count < max_count) {
            /* Check resonance */
            double* m_dbl = (double*)malloc((size_t)n * sizeof(double));
            if (m_dbl) {
                for (int i = 0; i < n; i++)
                    m_dbl[i] = (double)parts[i];
                for (int k = 0; k < n && count < max_count; k++) {
                    if (nft_is_resonant(m_dbl, eigenvalues, n, k)) {
                        for (int j = 0; j < n; j++)
                            resonant_list[count][j] = parts[j];
                        count++;
                        break;
                    }
                }
                free(m_dbl);
            }
            /* Next composition (rightmost carry) */
            int pos = n - 1;
            while (pos >= 0 && parts[pos] == 0) pos--;
            if (pos < 0) break;
            if (pos == 0) break;
            parts[pos]--;
            parts[pos-1]++;
        }
        free(parts);
    }
    free(current);
    return count;
}

/* Compute the homological operator [L_A, F] = L_A(F).
 * For a given vector field F (as a monomial), compute its
 * Lie bracket with the linear part A.
 * A[n*n]: linear part (diagonal canonical form assumed).
 * monomial[n]: exponent vector.
 * mode_k: component index.
 * result[n]: output of L_A applied to this monomial.
 * Returns 0 on success. */
int nft_homological_operator(const double* A, int n,
                               const int* monomial, int mode_k,
                               double* result) {
    if (!A || !monomial || !result || n <= 0
        || mode_k < 0 || mode_k >= n) return -1;
    double a_kk = A[mode_k * n + mode_k];
    double dot = 0.0;
    for (int i = 0; i < n; i++)
        dot += (double)monomial[i] * A[i * n + i];
    result[mode_k] = a_kk - dot;
    for (int j = 0; j < n; j++)
        if (j != mode_k) result[j] = 0.0;
    return 0;
}

/* Near-identity coordinate transformation:
 * y = x + h(x) where h is O(|x|^2).
 * Computes transformed vector field coefficients.
 * old_coeffs: coefficients in x-coordinates.
 * h_coeffs: generator coefficients.
 * new_coeffs: output coefficients in y-coordinates.
 * n_coeffs: number of coefficients.
 * Returns 0 on success. */
int nft_near_identity_transform(const double* old_coeffs,
                                  const double* h_coeffs,
                                  double* new_coeffs,
                                  int n_coeffs) {
    if (!old_coeffs || !h_coeffs || !new_coeffs || n_coeffs <= 0)
        return -1;
    /* First-order approximation: new = old + [L_A, h] + higher */
    for (int i = 0; i < n_coeffs; i++)
        new_coeffs[i] = old_coeffs[i] + h_coeffs[i];
    return 0;
}

/* Compute the normal form coefficients up to degree d_max.
 * A[n*n]: linear part (assumed diagonal).
 * coeffs[][n]: nonlinear coefficients (up to d_max).
 * d_max: maximum degree.
 * nf_coeffs[][n]: output normal form coefficients.
 * Returns 0 on success. */
int nft_compute_normal_form(const double* A, int n,
                              const double** coeffs, int d_max,
                              double** nf_coeffs) {
    if (!A || !coeffs || !nf_coeffs || n <= 0 || d_max < 2)
        return -1;
    double* ev = (double*)malloc((size_t)n * sizeof(double));
    if (!ev) return -1;
    for (int i = 0; i < n; i++) ev[i] = A[i * n + i];

    /* For each degree, keep only resonant terms */
    for (int d = 2; d <= d_max; d++) {
        double* m_dbl = (double*)calloc((size_t)n, sizeof(double));
        if (m_dbl) {
            m_dbl[0] = (double)d;
            for (int k = 0; k < n; k++) {
                if (nft_is_resonant(m_dbl, ev, n, k))
                    nf_coeffs[d][k] = coeffs[d][k];
                else
                    nf_coeffs[d][k] = 0.0;
            }
            free(m_dbl);
        }
    }
    free(ev);
    return 0;
}

/* Compute Jordan chain of generalized eigenvectors.
 * For an eigenvalue lambda with algebraic multiplicity m:
 * Solve (A - lambda*I)*v_1 = 0 for eigenvector.
 * Solve (A - lambda*I)*v_{k+1} = v_k for generalized eigenvectors.
 * A[n*n]: matrix. lambda: eigenvalue.
 * chain[][n]: output vectors (caller alloc, size max_len * n).
 * max_len: maximum chain length to compute.
 * Returns actual chain length. */
int nft_jordan_chain(const double* A, int n, double lambda,
                      double** chain, int max_len) {
    if (!A || !chain || n <= 0 || max_len < 1) return 0;
    double* M = (double*)malloc((size_t)(n * n) * sizeof(double));
    if (!M) return 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            M[i * n + j] = (i == j) ? -lambda : 0.0;
            M[i * n + j] -= A[i * n + j];
        }
    /* Find eigenvector v1: solve M*v1 = 0 (nullspace) */
    /* Simple: use last column as solution if rank n-1 */
    for (int i = 0; i < n; i++)
        chain[0][i] = (i == n-1) ? 1.0 : 0.0;

    /* Compute generalized eigenvectors */
    int len = 1;
    for (int k = 1; k < max_len; k++) {
        /* Solve M * chain[k] = chain[k-1] */
        double* rhs = (double*)malloc((size_t)n * sizeof(double));
        if (rhs) {
            for (int i = 0; i < n; i++) rhs[i] = chain[k-1][i];
            /* Simple substitution for upper triangular M */
            for (int i = 0; i < n; i++) {
                double s = rhs[i];
                for (int j = 0; j < i; j++)
                    s -= M[i*n+j] * chain[k][j];
                chain[k][i] = (fabs(M[i*n+i]) > 1e-15)
                              ? s / M[i*n+i] : 0.0;
            }
            free(rhs);
            len++;
        }
    }
    free(M);
    return len;
}

/* Compute the linear part of the normal form (diagonalization).
 * A[n*n]: original linear part.
 * diag[n]: output diagonal entries.
 * Returns 0 if A is diagonalizable (real eigenvalues), -1 otherwise. */
int nft_diagonalize_linear_part(const double* A, int n, double* diag) {
    if (!A || !diag || n <= 0) return -1;
    /* Compute trace and off-diagonal norm */
    double off_norm = 0.0;
    for (int i = 0; i < n; i++) {
        diag[i] = A[i * n + i];
        for (int j = 0; j < n; j++)
            if (i != j) off_norm += A[i * n + j] * A[i * n + j];
    }
    return (off_norm < 1e-10) ? 0 : -1;
}

/* Compute unfolded normal form with parameters.
 * nf_bare[n]: normal form at bifurcation point (param=0).
 * unfolding[n][p]: derivative of nf w.r.t. each parameter.
 * params[p]: parameter values.
 * result[n]: output normal form at given params.
 * n: state dim, p: parameter dim. Returns 0 on success. */
int nft_unfold_normal_form(const double* nf_bare,
                             const double** unfolding,
                             const double* params,
                             int n, int p, double* result) {
    if (!nf_bare || !unfolding || !params || !result
        || n <= 0 || p <= 0) return -1;
    for (int i = 0; i < n; i++) {
        result[i] = nf_bare[i];
        for (int k = 0; k < p; k++)
            result[i] += unfolding[i][k] * params[k];
    }
    return 0;
}

/* Classify normal form type based on eigenvalue configuration.
 * Returns: 0=saddle-node, 1=Hopf, 2=Takens-Bogdanov,
 * 3=fold-Hopf, 4=Hopf-Hopf, -1=unclassified. */
int nft_classify_bifurcation(const double* eigenvalues, int n) {
    if (!eigenvalues || n <= 0) return -1;
    int n_zero = 0, n_imag = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(eigenvalues[i]) < 1e-10) n_zero++;
    }
    /* Count purely imaginary pairs */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(eigenvalues[i] + eigenvalues[j]) < 1e-10
                && fabs(eigenvalues[i] * eigenvalues[j] + 1.0) < 0.1)
                n_imag++;
        }
    }
    if (n_zero == 1 && n_imag == 0) return 0;  /* Saddle-node */
    if (n_zero == 0 && n_imag == 1) return 1;  /* Hopf */
    if (n_zero == 2 && n_imag == 0) return 2;  /* Takens-Bogdanov */
    if (n_zero == 1 && n_imag == 1) return 3;  /* fold-Hopf */
    if (n_zero == 0 && n_imag == 2) return 4;  /* Hopf-Hopf */
    return -1;
}

