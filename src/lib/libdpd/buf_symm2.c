#include <stdio.h>
#include "dpd.h"

/* dpd_buf_symm2(): Symmetrizes two two-electron dpd buffers by
** taking, I'(pq,rs) = 1/2 [I1(pq,rs) + I2(pq,rs)] (note the
** indices!).  Users should keep in mind that the second buffer will
** be overwritten when this function is called.  Also note that this
** routine will NOT check to see if the row and column dimensions of
** the input buffers are identical, which is necessary for this to
** work.
**
** Arguments:
**   struct dpdbuf *Buf1: A pointer to the dpdbuf to be symmetrized.
**   struct dpdbuf *Buf2: A pointer to the dpdbuf to be symmetrized.  */

int dpd_buf_symm2(struct dpdbuf *Buf1, struct dpdbuf *Buf2)
{
  int h, row, col;
  double value;

  for(h=0; h < Buf1->params->nirreps; h++) {
      dpd_buf_mat_irrep_init(Buf1, h);
      dpd_buf_mat_irrep_rd(Buf1, h, 0, (FILE *) NULL);

      dpd_buf_mat_irrep_init(Buf2, h);
      dpd_buf_mat_irrep_rd(Buf2, h, 0, (FILE *) NULL);

      for(row=0; row < Buf1->params->rowtot[h]; row++)
	  for(col=0; col < Buf1->params->coltot[h]; col++) {
	      value = 0.5*(Buf1->matrix[h][row][col]+Buf2->matrix[h][col][row]);
	      Buf1->matrix[h][row][col] = value;
	    }

      dpd_buf_mat_irrep_wrt(Buf1, h, 0, (FILE *) NULL);
      dpd_buf_mat_irrep_close(Buf1, h);
      dpd_buf_mat_irrep_close(Buf2, h);
    }

  return 0;
}
