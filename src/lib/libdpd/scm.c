#include <stdio.h>
#include "dpd.h"

/* dpd_scm(): Multiplies every element of a two-electron dpdbuf by a scalar.
**
** Arguments:
**   struct dpdbuf *InBuf: A pointer to the dpdbuf.
**   double alpha: The scalar.
**   int print_flag: A booelan for the print routines.
**   FILE *outfile: A pointer to the formatted output stream.
*/

int dpd_scm(struct dpdbuf *InBuf, double alpha, int print_flag, FILE *outfile)
{
  int h, nirreps;
  int row, col;

  nirreps = InBuf->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(InBuf, h);
      dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);

      for(row=0; row < InBuf->params->rowtot[h]; row++) {
	  for(col=0; col < InBuf->params->coltot[h]; col++) {
	      InBuf->matrix[h][row][col] *= alpha;
	    }
	}

      dpd_buf_mat_irrep_wrt(InBuf, h, print_flag, outfile);
      dpd_buf_mat_irrep_close(InBuf, h);
    }

  return 0;
}
