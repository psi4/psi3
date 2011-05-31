/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

void transpert(const char *pert);

void roa(void)
{
  int h, h1, nirreps, row, col, dim;
  int Ga, Gi, i, a, ai, aa, ii;
  int j, k, l;
  int lwork, *ipiv, info;
  double *work;
  dpdbuf4 A, B;
  double **C;
  double **L, **R, **S, ***Q;
  double ***tensor_A, **tensor_alpha, **tensor_G;

  tensor_A = (double ***) malloc(3 * sizeof(double **));
  for(k=0; k < 3; k++) tensor_A[k] = block_matrix(3,3);

  tensor_alpha = block_matrix(3,3);
  tensor_G = block_matrix(3,3);

  Q = (double ***) malloc(3 * sizeof(double **));
  R = (double **) malloc(3 * sizeof(double *));
  L = (double **) malloc(3 * sizeof(double *));
  S = (double **) malloc(3 * sizeof(double *));

  dpd_buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_init(&B, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "B(AI,BJ)");

  for(h=0; h < moinfo.nirreps; h++) {

    dim = A.params->rowtot[h];

    if(dim) {

      dpd_buf4_mat_irrep_init(&A, h);
      dpd_buf4_mat_irrep_rd(&A, h);
      dpd_buf4_mat_irrep_init(&B, h);
      dpd_buf4_mat_irrep_rd(&B, h);

      C = block_matrix(2*dim, 2*dim);

      for(row=0; row < dim; row++) {
	for(col=0; col < dim; col++) {
	  C[row][col] = 2 * A.matrix[h][row][col];
	  C[row+dim][col+dim] = 2 * A.matrix[h][row][col];
	  C[row][col+dim] = -2 * B.matrix[h][row][col];
	  C[row+dim][col] = -2 * B.matrix[h][row][col];
	}

	C[row][row] += 2 * params.omega;
	C[row+dim][row+dim] -= 2 * params.omega;
      }

      ipiv = init_int_array(2*dim);
      lwork = 20 * 2*dim;
      work = init_array(lwork);
      info = C_DGETRF(2*dim, 2*dim, &(C[0][0]), 2*dim, ipiv);
      info = C_DGETRI(2*dim, &(C[0][0]), 2*dim, ipiv, work, lwork);
      if(info) {
	fprintf(outfile, "\n\tDGETRI failed. info = %d. Exiting.\n", info);
	exit(PSI_RETURN_FAILURE);
      }

      free(ipiv);
      free(work);

      transpert("Mu_X"); transpert("Mu_Y"); transpert("Mu_Z");
      for(k=0; k < 3; k++) {
        R[k] = init_array(2 * dim);
        for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
   	  Gi = h^Ga;
	  for(a=0; a < moinfo.virtpi[Ga]; a++) {
	    aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
	    for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
	      ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
	      R[k][ai] = 2 * moinfo.MU[k][aa][ii];
	      R[k][ai+dim] = -2 * moinfo.MU[k][aa][ii];
	    }
	  }
        }
      }

      transpert("L_X"); transpert("L_Y"); transpert("L_Z");
      for(k=0; k < 3; k++) {
        L[k] = init_array(2 * dim);
        for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
   	  Gi = h^Ga;
	  for(a=0; a < moinfo.virtpi[Ga]; a++) {
	    aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
	    for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
	      ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
	      L[k][ai] = 2 * moinfo.L[k][aa][ii];
	      L[k][ai+dim] = 2 * moinfo.L[k][aa][ii];
	    }
	  }
        }
      }

      transpert("Q_XX"); transpert("Q_XY"); transpert("Q_XZ");
      transpert("Q_YX"); transpert("Q_YY"); transpert("Q_YZ");
      transpert("Q_ZX"); transpert("Q_ZY"); transpert("Q_ZZ");

      Q[0] = block_matrix(3,2*dim);
      Q[1] = block_matrix(3,2*dim);
      Q[2] = block_matrix(3,2*dim);
      for(k=0; k < 3; k++)
        for(l=0; l < 3; l++) {

          for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
            Gi = h^Ga;
            for(a=0; a < moinfo.virtpi[Ga]; a++) {
              aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
              for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
                ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
                Q[k][l][ai] = 2 * moinfo.Q[k][l][aa][ii];
                Q[k][l][ai+dim] = -2 * moinfo.Q[k][l][aa][ii];
              }
            }
          }
        }

      for(k=0; k < 3; k++) {
        S[k] = init_array(2 * dim);
        C_DGEMV('t', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(R[k][0]), 1, 0, &(S[k][0]), 1);
      }

      for(i=0; i < 3; i++)
        for(j=0; j < 3; j++)
          for(k=0; k < 3; k++)
            for(row=0; row < 2*dim; row++)
	      tensor_A[i][j][k] += S[i][row] * Q[j][k][row];

      for(i=0; i < 3; i++)
        for(j=0; j < 3; j++)
          for(row=0; row < 2*dim; row++)
            tensor_alpha[i][j] += S[i][row] * R[j][row];

      for(i=0; i < 3; i++)
        for(j=0; j < 3; j++)
          for(row=0; row < 2*dim; row++)
            tensor_G[i][j] += S[i][row] * L[j][row];

      for(k=0; k < 3; k++) {
        free(S[k]);
        free(L[k]);
        free(R[k]);
        free_block(Q[k]);
      }

      free_block(C);

    } /* if(dim) */

    dpd_buf4_mat_irrep_close(&A, h);
    dpd_buf4_mat_irrep_close(&B, h);
  }

  fprintf(outfile, "\n\tHartree-Fock Alpha Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t----------------------------------------------\n");
  mat_print(tensor_alpha, 3, 3, outfile);

  fprintf(outfile, "\n\tHartree-Fock G' Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t-------------------------------------------\n");
  mat_print(tensor_G, 3, 3, outfile);

  fprintf(outfile, "\n\tHartree-Fock A Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t-------------------------------------------\n");
  for(k=0; k < 3; k++) {
    mat_print(tensor_A[k], 3, 3, outfile);
    free_block(tensor_A[k]);
  }
  free_block(tensor_alpha);
  free_block(tensor_G);

  dpd_buf4_close(&A);
  dpd_buf4_close(&B);
}

}} // namespace psi::response
