#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void diag_A_UHF(void)
{
  int h, i, dim, dim_A, dim_B;
  int nroot, root;
  int ai, bj, ck, c, k, C, K, Ksym;
  double **A, *eps, **v;
  char lbl[32];
  dpdbuf4 A_AA, A_BB, A_AB;
  dpdfile2 B;

  dpd_buf4_init(&A_AA, PSIF_MO_HESS, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
  dpd_buf4_init(&A_BB, PSIF_MO_HESS, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
  dpd_buf4_init(&A_AB, PSIF_MO_HESS, 0, 21, 31, 21, 31, 0, "A(AI,bj)");
  for(h=0; h < moinfo.nirreps; h++) {
    dim_A = A_AA.params->rowtot[h];
    dim_B = A_BB.params->rowtot[h];

    dim = dim_A + dim_B;
    moinfo.rank[h] = dim;

    A = block_matrix(dim, dim);
    eps = init_array(dim);
    v = block_matrix(dim, dim);

    dpd_buf4_mat_irrep_init(&A_AA, h);
    dpd_buf4_mat_irrep_rd(&A_AA, h);
    for(ai=0; ai < dim_A; ai++)
      for(bj=0; bj < dim_A; bj++)
	A[ai][bj] = A_AA.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_AA, h);

    dpd_buf4_mat_irrep_init(&A_BB, h);
    dpd_buf4_mat_irrep_rd(&A_BB, h);
    for(ai=0; ai < dim_B; ai++)
      for(bj=0; bj < dim_B; bj++)
	A[ai+dim_A][bj+dim_A] = A_BB.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_BB, h);

    dpd_buf4_mat_irrep_init(&A_AB, h);
    dpd_buf4_mat_irrep_rd(&A_AB, h);
    for(ai=0; ai < dim_A; ai++)
      for(bj=0; bj < dim_B; bj++)
	A[ai][bj+dim_A] = A[bj+dim_A][ai] = A_AB.matrix[h][ai][bj];
    dpd_buf4_mat_irrep_close(&A_AB, h);

    sq_rsp(dim, dim, A, eps, 1, v, 1e-12);

    for(i=0; i < MIN0(moinfo.rank[h],5); i++)
      moinfo.A_evals[h][i] = eps[i];

    free(eps);
    free_block(v);
    free_block(A);
  }

  dpd_buf4_close(&A_AA);
  dpd_buf4_close(&A_BB);
  dpd_buf4_close(&A_AB);

}
