#include <stdio.h>
#include <qt.h>
#include "dpd.h"

double dpd_oe_dot(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		  int print_flag, FILE *outfile)
{
  int h, nirreps;
  double dot;

  nirreps = FileA->params->nirreps;

  dot = 0.0;

  for(h=0; h < nirreps; h++) {

      dpd_oe_file_mat_irrep_init(FileA, h);
      dpd_oe_file_mat_irrep_init(FileB, h);
      dpd_oe_file_mat_irrep_rd(FileA, h, print_flag, outfile);
      dpd_oe_file_mat_irrep_rd(FileB, h, print_flag, outfile);

      dot += dot_block(FileA->matrix[h], FileB->matrix[h],
		       FileA->params->rowtot[h], FileA->params->coltot[h], 1.0);
      
      dpd_oe_file_mat_irrep_close(FileA, h);
      dpd_oe_file_mat_irrep_close(FileB, h);
    }

  return dot;

}
