#include <stdio.h>
#include <math.h>
#include <qt.h>
#include "dpd.h"

/* dpd_contract121(): Contracts a two-electron dpd with a one-electron
** dpd to give another one-electron dpd.  Both indices of the
** one-electron factor must be summed, and both must be ket indices in
** the two-electron buffer.
**
** Arguments:
**   struct dpdbuf *X: A pointer to the two-electron buffer.
**   struct oe_dpdfile *Y: A pointer to the one-electron factor file.
**   struct oe_dpdfile *Z: A pointer to the one-electron target file.
**   int trans_Y: A boolean to indicate whether the indices in Y are
**                transposed relative to those in the ket of X.
**   int trans_Z: A boolean to indicate whether the indices in Z are
**                transposed relative to those in the bra of X.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_contract121(struct dpdbuf *X, struct oe_dpdfile *Y, struct
		    oe_dpdfile *Z, int trans_Y, int trans_Z, double
		    alpha, double beta, int print_flag, FILE *outfile)
{
  int nirreps;
  int row,p,q,r,s, psym, qsym, Gr, Gs, P, Q, R, S, col;
  double **TMP;
  double value;
#ifdef DPD_DEBUG
  int *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  dpd_oe_file_mat_init(Y);
  dpd_oe_file_mat_rd(Y, print_flag, outfile);
  dpd_oe_file_mat_init(Z);
  if(fabs(beta) > 0.0) dpd_oe_file_mat_rd(Z, print_flag, outfile);

#ifdef DPD_DEBUG
  if(trans_Z) { zrow = Z->params->coltot; zcol = Z->params->rowtot; }
  else { zrow = Z->params->rowtot; zcol = Z->params->coltot; }

  if(trans_Y) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

  if((zrow != X->params->ppi) || (zcol != X->params->qpi) ||
     (yrow != X->params->rpi) || (ycol != X->params->spi)) {
      fprintf(outfile, "** Alignment error in contract121 **\n");
      dpd_error("dpd_contract121", outfile);
    }
#endif

  dpd_buf_mat_irrep_init(X, 0);
  dpd_buf_mat_irrep_rd(X, 0, print_flag, outfile);

  for(row=0; row < X->params->rowtot[0]; row++) {
      p = X->params->roworb[0][row][0];
      psym = X->params->psym[p];
      P = p - X->params->poff[psym];
      q = X->params->roworb[0][row][1];
      qsym = X->params->qsym[q];
      Q = q - X->params->qoff[qsym];

      value = 0.0;

      for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr;

	  if(X->params->spi[Gs] && X->params->rpi[Gr]) {
	      if(trans_Y) {
		  TMP = block_matrix(X->params->spi[Gs],X->params->rpi[Gr]);
		}
	      else {
		  TMP = block_matrix(X->params->rpi[Gr],X->params->spi[Gs]);
		}
	    }

	  for(r=0; r < X->params->rpi[Gr]; r++) {
	      R = X->params->roff[Gr] + r;
	      for(s=0; s < X->params->spi[Gs]; s++) {
		  S = X->params->soff[Gs] + s;

		  col = X->params->colidx[R][S];

		  if(trans_Y) 
		      TMP[s][r] = X->matrix[0][row][col];
		  else
		      TMP[r][s] = X->matrix[0][row][col];
		}
	    }

	  if(trans_Y) 
	      value += dot_block(TMP, Y->matrix[Gs], X->params->spi[Gs],
				 X->params->rpi[Gr], alpha);
	  else
	      value += dot_block(TMP, Y->matrix[Gr], X->params->rpi[Gr],
				 X->params->spi[Gs], alpha);

	  if(X->params->rpi[Gr] && X->params->spi[Gs])
	      free_block(TMP);
	    }

      if(trans_Z) 
	  Z->matrix[qsym][Q][P] = beta*Z->matrix[qsym][Q][P] + value;
      else 
	  Z->matrix[psym][P][Q] = beta*Z->matrix[psym][P][Q] + value;
    }

  dpd_buf_mat_irrep_close(X, 0);

  dpd_oe_file_mat_close(Y);
  dpd_oe_file_mat_wrt(Z, print_flag, outfile);
  dpd_oe_file_mat_close(Z);

  return 0;
}
