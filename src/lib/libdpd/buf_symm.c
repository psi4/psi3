#include <stdio.h>
#include "dpd.h"

/* dpd_buf_symm(): Symmetrizes a two-electron dpd buffer by taking,
** I'(pq,rs) = 1/2 [I(pq,rs) + I(rs,pq)].  Users should keep in mind
** that the original buffer will be overwritten when this function is
** called.  Also note that this routine will NOT check to see if the
** row and column dimensions of the input buffer are identical, which
** is necessary for this to work.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the dpdbuf to be symmetrized.
*/

int dpd_buf_symm(struct dpdbuf *Buf)
{
  int h, row, col;
  double value;

  for(h=0; h < Buf->params->nirreps; h++) {
      dpd_buf_mat_irrep_init(Buf, h);
      dpd_buf_mat_irrep_rd(Buf, h, 0, (FILE *) NULL);

      for(row=0; row < Buf->params->rowtot[h]; row++)
	  for(col=0; col < Buf->params->coltot[h]; col++) {
	      value = 0.5*(Buf->matrix[h][row][col]+Buf->matrix[h][col][row]);
	      Buf->matrix[h][row][col] = Buf->matrix[h][col][row] = value;
	    }

      dpd_buf_mat_irrep_wrt(Buf, h, 0, (FILE *) NULL);
      dpd_buf_mat_irrep_close(Buf, h);
    }

  return 0;
}
