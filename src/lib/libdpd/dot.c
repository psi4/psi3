#include <stdio.h>
#include <qt.h>
#include "dpd.h"

double dpd_dot(struct dpdbuf *BufA, struct dpdbuf *BufB, 
	       int print_flag, FILE *outfile)
{
  int h, nirreps;
  double dot;

  nirreps = BufA->params->nirreps;

  dot = 0.0;

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(BufA, h);
      dpd_buf_mat_irrep_init(BufB, h);
      dpd_buf_mat_irrep_rd(BufA, h, print_flag, outfile);
      dpd_buf_mat_irrep_rd(BufB, h, print_flag, outfile);

      dot += dot_block(BufA->matrix[h], BufB->matrix[h],
		       BufA->params->rowtot[h], BufA->params->coltot[h], 1.0);

      dpd_buf_mat_irrep_close(BufA, h);
      dpd_buf_mat_irrep_close(BufB, h);
    }

  return dot;
}
      
