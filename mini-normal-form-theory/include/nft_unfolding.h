#ifndef NFT_UNFOLDING_H
#define NFT_UNFOLDING_H
#include "nft_types.h"
Unfolding* nft_unfolding_create(int n_params);
void nft_unfolding_free(Unfolding* uf);
int nft_compute_versal_unfolding(const NFTMatrix* nf, Unfolding* uf);
bool nft_unfolding_is_versal(const Unfolding* uf);
double nft_unfolding_codimension(const NFTMatrix* nf);
BifurcationType nft_classify_bifurcation(const NormalFormResult* nf, double param);
const char* nft_bifurcation_name(BifurcationType bt);
int nft_unfolding_evaluate(const Unfolding* uf, const double* x, const double* mu, double* dxdt);
void nft_unfolding_print(const Unfolding* uf);
int nft_tangent_space_dimension(const NFTMatrix* A, const NormalFormResult* nf);
int nft_miniversal_unfolding(const NFTMatrix* nf, Unfolding* uf);
void nft_unfolding_jacobian(const Unfolding* uf, const double* x, const double* mu, NFTMatrix* J);
int nft_bifurcation_diagram(const Unfolding* uf, double mu_min, double mu_max, int n_steps, double* eq, int* n_eq);
int nft_detect_bifurcation_point(const Unfolding* uf, double mu_start, double mu_end, double* mu_crit);
BifurcationType nft_classify_from_unfolding(const Unfolding* uf, double mu_range);
int nft_lyapunov_schmidt_bifurcation(const NFTMatrix* A, double* g, double* dg, int* dim_null);
#endif
