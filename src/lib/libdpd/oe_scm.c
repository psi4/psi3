#include <stdio.h>
#include "dpd.h"

int dpd_oe_scm(struct oe_dpdfile *InFile, double alpha, int print_flag,
	       FILE *outfile)
{
  int h, nirreps;
  int row, col;

  nirreps = InFile->params->nirreps;

  dpd_oe_file_mat_init(InFile);
  dpd_oe_file_mat_rd(InFile, print_flag, outfile);

  for(h=0; h < nirreps; h++) {

      for(row=0; row < InFile->params->rowtot[h]; row++) {
	  for(col=0; col < InFile->params->coltot[h]; col++) {
	      InFile->matrix[h][row][col] *= alpha;
	    }
	}
    }

  dpd_oe_file_mat_wrt(InFile, print_flag, outfile);
  dpd_oe_file_mat_close(InFile);

  return 0;
}
