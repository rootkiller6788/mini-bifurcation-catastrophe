#ifndef NFT_CENTER_H
#define NFT_CENTER_H
#include "nft_types.h"
CenterManifold* nft_center_create(int dim_c, int dim_s, int order);
void nft_center_free(CenterManifold* cm);
int nft_center_compute(CenterManifold* cm, const NFTMatrix* A, const NFTMatrix* f_nl);
int nft_center_compute_to_order(CenterManifold* cm, const NFTMatrix* A, const NFTMatrix* f_nl, int target_order);
void nft_center_evaluate(const CenterManifold* cm, const double* x_c, double* x_s);
int nft_center_reduced_system(const CenterManifold* cm, const NFTMatrix* A, NFTMatrix** A_red);
double nft_center_approximation_error(const CenterManifold* cm, const NFTMatrix* A);
void nft_center_print(const CenterManifold* cm);
double nft_center_eval_taylor(const CenterManifold* cm, const double* x_c, int component, int order);
bool nft_center_verify_tangency(const CenterManifold* cm);
void nft_center_reduced_dynamics(const CenterManifold* cm, const NFTMatrix* A, const double* x_c, double* dx_c);
int nft_split_eigenvalues(const double* ev, int n, int* dim_c, int* dim_s, int* dim_u, double tol);
int nft_transform_to_center_stable(const NFTMatrix* A, NFTMatrix** T, int* dim_c, int* dim_s, int* dim_u);
int nft_spectral_decomposition(const NFTMatrix* A, NFTMatrix** Ac, NFTMatrix** As, NFTMatrix** Au, double tol);
int nft_lyapunov_schmidt_dimension(const NFTMatrix* A, double tol);
int nft_center_suspended(CenterManifold* cm, const NFTMatrix* A, double param, int* new_dim_c);
int nft_center_normal_form(const CenterManifold* cm, NFTMatrix** nf_red);
int nft_center_stability_type(const CenterManifold* cm, const NFTMatrix* A);
void nft_print_spectral_splitting(const NFTMatrix* A, double tol);
#endif
