#include <stdio.h>
#include "dpd.h"

/* dpd_axpy(): Evaluates the standard operation a * X + Y -> Y for dpd
** two-electron buffers.
**
** Arguments:
**   struct dpdbuf *BufX: A pointer to the leftmost dpdbuf.
**   struct dpdbuf *BufY: A pointer to the rightmost (and target)
**                        dpdbuf.
**   double alpha: The scalar prefactor in the multiplication.
**   int print_flag: A booelan for the print routines.
**   FILE *outfile: A pointer to the formatted output stream.
*/

int dpd_axpy(struct dpdbuf *BufX, struct dpdbuf *BufY, double alpha,
	     int print_flag, FILE *outfile)
{
  int h, nirreps;
  int row, col;

  nirreps = BufX->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(BufX, h);
      dpd_buf_mat_irrep_rd(BufX, h, print_flag, outfile);

      dpd_buf_mat_irrep_init(BufY, h);
      dpd_buf_mat_irrep_rd(BufY, h, print_flag, outfile);

      for(row=0; row < BufX->params->rowtot[h]; row++) {
	  for(col=0; col < BufX->params->coltot[h]; col++) {
	      BufY->matrix[h][row][col] += alpha*BufX->matrix[h][row][col]; 
	    }
	}

      dpd_buf_mat_irrep_wrt(BufY, h, print_flag, outfile);

      dpd_buf_mat_irrep_close(BufX, h);
      dpd_buf_mat_irrep_close(BufY, h);
    }

  return 0;
}
