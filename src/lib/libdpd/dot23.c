#include <stdio.h>
#include <qt.h>
#include "dpd.h"

/* dpd_dot23(): Contracts a one-electron dpd with a two-electron dpd
** where both indices in the former and indices two and three in the
** latter are summed.
**
** Arguments:
**   struct oe_dpdfile *T: A pointer to the one-electron factor.
**   struct dpdbuf *I: A pointer to the two-electron factor.
**   struct oe_dpdfile *Z: A pointer to the one-electron target.
**   int transt: A boolean indicating whether the T-factor should be
**               transposed.
**   int transz: A boolean indicating whether the Z-product should be
**               transposed.
**   double alpha: A prefactor for the product alpha * T * I.
**   double beta: A prefactor for the target beta * Z.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
**
** Example contractions:
**    beta * Z(p,s) = alpha * T(q,r) * I(pq,rs)
**           (transt = 0; transz =0;)
**    beta * Z(p,s) = alpha * T(r,q) * I(pq,rs)
**           (transt = 1; transz =0;)
**    beta * Z(s,p) = alpha * T(q,r) * I(pq,rs)
**           (transt = 0; transz =1;)
**    beta * Z(s,p) = alpha * T(r,q) * I(pq,rs)
**           (transt = 1; transz =1;)
*/
   

int dpd_dot23(struct oe_dpdfile *T, struct dpdbuf *I, struct oe_dpdfile *Z,
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
	  Gs = Gp;  Gq = Gr = h^Gp;

	  /* Allocate space for the X buffer */
	  if(T->params->ppi[Gq] && T->params->qpi[Gr])
	      X = block_matrix(T->params->ppi[Gq],T->params->qpi[Gr]);

	  /* Loop over orbitals of the target */
	  for(p=0; p < Z->params->ppi[Gp]; p++) {
	      P = Z->params->poff[Gp] + p;
	      for(s=0; s < Z->params->qpi[Gs]; s++) {
		  S = Z->params->qoff[Gs] + s;

		  /* Loop over orbitals of the two-index term */
		  for(q=0; q < T->params->ppi[Gq]; q++) {
		      Q = T->params->poff[Gq] + q;
		      for(r=0; r < T->params->qpi[Gr]; r++) {
			  R = T->params->qoff[Gr] + r;

			  /* Calculate row and column indices in I */
                          if(!transt && !transz) {
                              row = I->params->rowidx[P][Q];
                              col = I->params->colidx[R][S];
                            }
                          else if(transt && !transz) {
                              row = I->params->rowidx[P][R];
                              col = I->params->colidx[Q][S];
                            }
                          else if(!transt && transz) {
                              row = I->params->rowidx[S][Q];
                              col = I->params->colidx[R][P];
                            }
                          else if(transt && transz) {
                              row = I->params->rowidx[S][R];
                              col = I->params->colidx[Q][P];
                            }

			  /* Build the X buffer */
			  X[q][r] = I->matrix[h][row][col]; 

			}
		    }

		  value = dot_block(T->matrix[Gq], X, T->params->ppi[Gq],
				    T->params->qpi[Gr], alpha); 
		  
		  Z->matrix[Gp][p][s] += value;
		}
	    }
	  if(T->params->ppi[Gq] && T->params->qpi[Gr])
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
