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

/* polar(): Compute the RHF dipole-polarizability tensor.

 In build_A_RHF(), the A(ai,bj) matrix was defined as:

  A(ai,bj) = delta_ij f_ab - delta_ab f_ij + 2 <ij|ab> - <ia|jb>

 In build_B_RHF(), the B(ai,bj) matrix was defined as:

  B(ai,bj) = 2<ij|ab> - <ij|ba>

 Here we modify these in agreement with Jorgensen and Simons (1981), p.120

  A(ai,bj) = 2[ delta_ij f_ab - delta_ab f_ij + 2 <ij|ab> - <ia|jb>]

  and

  B(ai,bj) = 2[<ij|ba> - 2<ij|ab>]

*/

namespace psi { namespace response {

void transpert(const char *pert);

void polar(void)
{
  int h, h1, nirreps, row, col, dim;
  int Ga, Gi, i, a, ai, aa, ii;
  int lwork, *ipiv, info;
  double *work;
  dpdbuf4 A, B;
  double **C;
  double *Rx, *Ry, *Rz, *Sx, *Sy, *Sz;
  double **alpha;
  double xnorm, ynorm, znorm;

  alpha = block_matrix(3,3); /* polarizability tensor */
  transpert("Mu_X");
  transpert("Mu_Y");
  transpert("Mu_Z");

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
	  C[row+dim][col+dim] = 2 *A.matrix[h][row][col];
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
	    Rx[ai+dim] = - 2 * moinfo.MU[0][aa][ii];
	    Ry[ai] = 2 * moinfo.MU[1][aa][ii];
	    Ry[ai+dim] = - 2 * moinfo.MU[1][aa][ii];
	    Rz[ai] = 2 * moinfo.MU[2][aa][ii];
	    Rz[ai+dim] = - 2 * moinfo.MU[2][aa][ii];
	  }
	}
      }

      Sx = init_array(2 * dim);
      Sy = init_array(2 * dim);
      Sz = init_array(2 * dim);

      C_DGEMV('n', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Rx[0]), 1, 0, &(Sx[0]), 1);
      C_DGEMV('n', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Ry[0]), 1, 0, &(Sy[0]), 1);
      C_DGEMV('n', 2*dim, 2*dim, 1, &(C[0][0]), 2*dim, &(Rz[0]), 1, 0, &(Sz[0]), 1);

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

	alpha[0][0] += Rx[row] * Sx[row];
	alpha[1][0] += Ry[row] * Sx[row];
	alpha[2][0] += Rz[row] * Sx[row];

	alpha[1][1] += Ry[row] * Sy[row];
	alpha[2][1] += Rz[row] * Sy[row];

	alpha[2][2] += Rz[row] * Sz[row];
      }

      free(Sx); free(Sy); free(Sz);
      free(Rx); free(Ry); free(Rz);

      free_block(C);

    } /* if(dim) */

    dpd_buf4_mat_irrep_close(&A, h);
    dpd_buf4_mat_irrep_close(&B, h);
  }

  /* fill out the tensor */
  alpha[0][1] = alpha[1][0];
  alpha[0][2] = alpha[2][0];
  alpha[1][2] = alpha[2][1];

  fprintf(outfile, "\n\tHartree-Fock Electric Polarizability Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t---------------------------------------------------------------\n");
  mat_print(alpha, 3, 3, outfile);
  free_block(alpha);

  dpd_buf4_close(&A);
  dpd_buf4_close(&B);
}

}} // namespace psi::response
