#include <stdio.h>
#include "dpd.h"

/* dpd_chksum(): Evaluates the sum of the squares of the elements of a
** given dpd two-electron buffer.
**
** Arguments:
**   struct dpdbuf *BufX: A pointer to the leftmost dpdbuf.
*/

double dpd_chksum(struct dpdbuf *BufX)
{
  int h, nirreps;
  int row, col;
  double alpha=0.0;

  nirreps = BufX->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(BufX, h);
      dpd_buf_mat_irrep_rd(BufX, h, 0, NULL);

      for(row=0; row < BufX->params->rowtot[h]; row++)
	  for(col=0; col < BufX->params->coltot[h]; col++)
              alpha += BufX->matrix[h][row][col] * BufX->matrix[h][row][col];

      dpd_buf_mat_irrep_close(BufX, h);
    }

  return alpha;
}
