#include <stdio.h>
#include "dpd.h"
#include <qt.h>

/* dpd_buf4_axpy(): Evaluates the standard operation a * X + Y -> Y for dpd
** four-index buffers.
**
** Arguments:
**   dpdbuf4 *BufX: A pointer to the leftmost dpdbuf4.
**   dpdbuf4 *BufY: A pointer to the rightmost (and target)
**                        dpdbuf4.
**   double alpha: The scalar prefactor in the multiplication.
*/

int dpd_buf4_axpy(dpdbuf4 *BufX, dpdbuf4 *BufY, double alpha)
{
  int h, nirreps, my_irrep;
  int row, col, length;
  double *X, *Y;

  nirreps = BufX->params->nirreps;
  my_irrep = BufX->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_axpy");
#endif

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(BufX, h);
      dpd_buf4_mat_irrep_rd(BufX, h);

      dpd_buf4_mat_irrep_init(BufY, h);
      dpd_buf4_mat_irrep_rd(BufY, h);

      /* I need to replace this with a BLAS1 call */
/*
      for(row=0; row < BufX->params->rowtot[h]; row++) {
	  for(col=0; col < BufX->params->coltot[h^my_irrep]; col++) {
	      BufY->matrix[h][row][col] += alpha*BufX->matrix[h][row][col]; 
	    }
	}
*/

      length = BufX->params->rowtot[h] * BufX->params->coltot[h];
      if(length) {
          X = &(BufX->matrix[h][0][0]);
          Y = &(BufY->matrix[h][0][0]);
          C_DAXPY(length, alpha, X, 1, Y, 1);
       }

      dpd_buf4_mat_irrep_wrt(BufY, h);

      dpd_buf4_mat_irrep_close(BufX, h);
      dpd_buf4_mat_irrep_close(BufY, h);
    }

#ifdef DPD_TIMER
  timer_off("buf4_axpy");
#endif

  return 0;
}
