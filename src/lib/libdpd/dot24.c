#include <stdio.h>
#include <qt.h>
#include "dpd.h"

int dpd_dot24(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
	      int transt, int transz, double alpha, double beta,
	      int print_flag, FILE *outfile)
{
  int h, Gp, Gq, Gr, Gs;
  int p, q, r, s;
  int P, Q, R, S;
  int row, col;
  int nirreps;
  double **X;
  double value;

  nirreps = T->params->nirreps;

  /* Get the two-index quantities from disk */
  dpd_oe_file_mat_init(T);
  dpd_oe_file_mat_rd(T, print_flag, outfile);
  dpd_oe_scm(Z, beta, print_flag, outfile);
  dpd_oe_file_mat_init(Z);
  dpd_oe_file_mat_rd(Z, print_flag, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(I, h);
      dpd_buf_mat_irrep_rd(I, h, print_flag, outfile);

      /* Loop over irreps of the target */
      for(Gp=0; Gp < nirreps; Gp++) {
	  Gr = Gp;  Gq = Gs = h^Gp;

	  /* Allocate space for the X buffer */
	  if(T->params->ppi[Gq] && T->params->qpi[Gs])
	      X = block_matrix(T->params->ppi[Gq],T->params->qpi[Gs]);

	  /* Loop over orbitals of the target */
	  for(p=0; p < Z->params->ppi[Gp]; p++) {
	      P = Z->params->poff[Gp] + p;
	      for(r=0; r < Z->params->qpi[Gr]; r++) {
		  R = Z->params->qoff[Gr] + r;

		  /* Loop over orbitals of the two-index term */
		  for(q=0; q < T->params->ppi[Gq]; q++) {
		      Q = T->params->poff[Gq] + q;
		      for(s=0; s < T->params->qpi[Gs]; s++) {
			  S = T->params->qoff[Gs] + s;

			  /* Calculate row and column indices in I */
                          if(!transt && !transz) {
                              row = I->params->rowidx[P][Q];
                              col = I->params->colidx[R][S];
                            }
                          else if(transt && !transz) {
                              row = I->params->rowidx[P][S];
                              col = I->params->colidx[R][Q];
                            }
                          else if(!transt && transz) {
                              row = I->params->rowidx[R][Q];
                              col = I->params->colidx[P][S];
                            }
                          else if(transt && transz) {
                              row = I->params->rowidx[R][S];
                              col = I->params->colidx[P][Q];
                            }

			  /* Build the X buffer */
			  X[q][s] = I->matrix[h][row][col]; 

			}
		    }

		  value = dot_block(T->matrix[Gq], X, T->params->ppi[Gq],
				    T->params->qpi[Gs], alpha); 
		  
		  Z->matrix[Gp][p][r] += value;
		}
	    }
	  if(T->params->ppi[Gq] && T->params->qpi[Gs])
	      free_block(X);
	}
      dpd_buf_mat_irrep_close(I, h);
    }

  /* Close the two-index quantities */
  dpd_oe_file_mat_close(T);
  dpd_oe_file_mat_wrt(Z, print_flag, outfile);
  dpd_oe_file_mat_close(Z);

  return 0;
}
