#include <stdio.h>
#include <qt.h>
#include "dpd.h"

double dpd_buf4_dot(dpdbuf4 *BufA, dpdbuf4 *BufB)
{
  int h, nirreps, my_irrep;
  double dot;

  nirreps = BufA->params->nirreps;
  my_irrep = BufA->file.my_irrep;

  dot = 0.0;

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(BufA, h);
      dpd_buf4_mat_irrep_init(BufB, h);
      dpd_buf4_mat_irrep_rd(BufA, h);
      dpd_buf4_mat_irrep_rd(BufB, h);

      dot += dot_block(BufA->matrix[h], BufB->matrix[h],
		       BufA->params->rowtot[h],
		       BufA->params->coltot[h^my_irrep], 1.0);

      dpd_buf4_mat_irrep_close(BufA, h);
      dpd_buf4_mat_irrep_close(BufB, h);
    }

  return dot;
}
      
