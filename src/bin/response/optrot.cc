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

void optrot(void)
{
  int h, h1, nirreps, row, col, dim;
  int Ga, Gi, i, a, ai, aa, ii;
  int lwork, *ipiv, info;
  double *work;
  dpdbuf4 A, B;
  double **C;
  double *Rx, *Ry, *Rz, *Sx, *Sy, *Sz, *Lx, *Ly, *Lz;
  double **alpha;
  double xnorm, ynorm, znorm;

  alpha = block_matrix(3,3); /* G' tensor */

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

      /* Set up the dipole vectors for this irrep */
      Rx = init_array(2 * dim);
      Ry = init_array(2 * dim);
      Rz = init_array(2 * dim);
      for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
	Gi = h^Ga;
	for(a=0; a < moinfo.virtpi[Ga]; a++) {
	  aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
	  for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
	    ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
	    Rx[ai] = 2 * moinfo.MU[0][aa][ii];
	    Rx[ai+dim] = -2 * moinfo.MU[0][aa][ii];
	    Ry[ai] = 2 * moinfo.MU[1][aa][ii];
	    Ry[ai+dim] = -2 * moinfo.MU[1][aa][ii];
	    Rz[ai] = 2 * moinfo.MU[2][aa][ii];
	    Rz[ai+dim] = -2 * moinfo.MU[2][aa][ii];
	  }
	}
      }

      transpert("L*_X"); transpert("L*_Y"); transpert("L*_Z");

      Lx = init_array(2 * dim);
      Ly = init_array(2 * dim);
      Lz = init_array(2 * dim);
      for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
        Gi = h^Ga;
        for(a=0; a < moinfo.virtpi[Ga]; a++) {
          aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
          for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
            ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
            Lx[ai] = 2 * moinfo.L[0][aa][ii];
            Lx[ai+dim] = -2 * moinfo.L[0][aa][ii];
            Ly[ai] = 2 * moinfo.L[1][aa][ii];
            Ly[ai+dim] = -2 * moinfo.L[1][aa][ii];
            Lz[ai] = 2 * moinfo.L[2][aa][ii];
            Lz[ai+dim] = -2 * moinfo.L[2][aa][ii];
          }
        }
      }

      Sx = init_array(2 * dim);
      Sy = init_array(2 * dim);
      Sz = init_array(2 * dim);

      C_DGEMV('t', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Rx[0]), 1, 0, &(Sx[0]), 1);
      C_DGEMV('t', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Ry[0]), 1, 0, &(Sy[0]), 1);
      C_DGEMV('t', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Rz[0]), 1, 0, &(Sz[0]), 1);

      xnorm = ynorm = znorm = 0.0;
      for(row=0; row < 2*dim; row++) {
       xnorm += Sx[0] * Sx[0];
       ynorm += Sy[0] * Sy[0];
       znorm += Sz[0] * Sz[0];
      }

      xnorm = sqrt(xnorm);
      ynorm = sqrt(ynorm);
      znorm = sqrt(znorm);
      fprintf(outfile, "XNorm = %20.10f\n", xnorm);
      fprintf(outfile, "YNorm = %20.10f\n", ynorm);
      fprintf(outfile, "ZNorm = %20.10f\n", znorm);

      for(row=0; row< 2*dim; row++) {

	alpha[0][0] += Lx[row] * Sx[row];
	alpha[1][0] += Ly[row] * Sx[row];
	alpha[0][1] += Lx[row] * Sy[row];
	alpha[0][2] += Lx[row] * Sz[row];
	alpha[2][0] += Lz[row] * Sx[row];
	alpha[1][1] += Ly[row] * Sy[row];
	alpha[1][2] += Ly[row] * Sz[row];
	alpha[2][1] += Lz[row] * Sy[row];
	alpha[2][2] += Lz[row] * Sz[row];
      }

      free(Sx); free(Sy); free(Sz);
      free(Lx); free(Ly); free(Lz);
      free(Rx); free(Ry); free(Rz);

      free_block(C);

    } /* if(dim) */

    dpd_buf4_mat_irrep_close(&A, h);
    dpd_buf4_mat_irrep_close(&B, h);
  }

  fprintf(outfile, "\n\tHartree-Fock G' Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t-------------------------------------------\n");
  mat_print(alpha, 3, 3, outfile);
  free_block(alpha);

  dpd_buf4_close(&A);
  dpd_buf4_close(&B);
}

}} // namespace psi::response
