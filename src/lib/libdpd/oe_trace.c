#include <stdio.h>
#include "dpd.h"

double dpd_oe_trace(struct oe_dpdfile *InFile, int print_flag, FILE *outfile)
{
  int h, nirreps;
  int row, col;
  double trace;

  nirreps = InFile->params->nirreps;

  dpd_oe_file_mat_init(InFile);
  dpd_oe_file_mat_rd(InFile, print_flag, outfile);

  trace = 0.0;
  for(h=0; h < nirreps; h++)
      for(row=0; row < InFile->params->rowtot[h]; row++)
	  trace += InFile->matrix[h][row][row];

  dpd_oe_file_mat_wrt(InFile, print_flag, outfile);
  dpd_oe_file_mat_close(InFile);

  return trace;
}
