#include <stdio.h>
#include "dpd.h"

/* dpd_oe_chksum(): Evaluates the sum of the squares of the elements of a
** given dpd one-electron buffer.
*/

double dpd_oe_chksum(struct oe_dpdfile *BufX)
{
  int h, nirreps;
  int row, col;
  double alpha=0.0;

  nirreps = BufX->params->nirreps;

  dpd_oe_file_mat_init(BufX);
  dpd_oe_file_mat_rd(BufX, 0, NULL);

  for(h=0; h < nirreps; h++) {

      for(row=0; row < BufX->params->rowtot[h]; row++)
	  for(col=0; col < BufX->params->coltot[h]; col++)
              alpha += BufX->matrix[h][row][col] * BufX->matrix[h][row][col];

    }

  dpd_oe_file_mat_close(BufX);

  return alpha;
}
