#ifndef NFT_TYPES_H
#define NFT_TYPES_H
#include <stdbool.h>
#include <stddef.h>

/* Normal Form Theory — Core Types
 * Poincare-Dulac normal forms, center manifold reduction,
 * versal unfoldings, homological equation.
 * Ref: Arnold (1983), Guckenheimer & Holmes (1983), Wiggins (2003) */

typedef struct { double* data; int rows; int cols; bool owns; } NFTMatrix;
typedef struct { double real; double imag; } NFTComplex;

/* Jordan block structure */
typedef enum { NFT_JORDAN_REAL=0, NFT_JORDAN_COMPLEX=1, NFT_JORDAN_SEMISIMPLE=2 } JordanType;
typedef struct { JordanType type; int size; double eigenvalue_real; double eigenvalue_imag; } JordanBlock;

/* Bifurcation classification */
typedef enum { BIF_NONE=0, BIF_SADDLE_NODE=1, BIF_HOPF=2, BIF_PITCHFORK=3, BIF_TRANSCRITICAL=4, BIF_HOPF_DEGENERATE=5, BIF_BOGDANOV_TAKENS=6, BIF_ZERO_HOPF=7 } BifurcationType;
typedef enum { UNFOLDING_MINIMAL=0, UNFOLDING_VERSAL=1, UNFOLDING_OVERFULL=2 } UnfoldingType;

/* Normal form specification */
typedef struct { int* orders; int n_vars; int n_params; double bifurcation_param; BifurcationType bif_type; } NormalFormSpec;

/* Resonance condition: lambda*m == lambda_j for some integer vector m with sum >= 2 */
typedef struct { int* m; int n; int sum; bool is_resonant; } ResonanceMonomial;

/* Homological equation: L_A(P) = Q where L_A is the Lie bracket */
typedef struct { NFTMatrix* A; NFTMatrix* Q; NFTMatrix* P; int n; bool solved; double residual; } HomologicalEquation;

/* Center manifold: W^c = { (x,y): y = h(x) } */
typedef struct { NFTMatrix* h_coeff; int dim_center; int dim_stable; int order; bool computed; } CenterManifold;

/* Versal unfolding: F(x,mu) = v(x) + sum_i alpha_i * e_i(x) */
typedef struct { NFTMatrix** basis; int n_basis; double* parameters; int n_params; bool is_versal; } Unfolding;

/* Normal form result */
typedef struct { NFTMatrix* normal_form; NFTMatrix* transformation; int* resonant_orders; int n_resonant; double truncation_error; int iterations; bool converged; } NormalFormResult;
#endif
