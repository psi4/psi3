#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_file_mat_irrep_row_zero(struct dpdfile *File, int irrep, int row)
{
  int coltot;
  
  coltot = File->params->coltot[irrep];

  if(coltot)
      zero_arr(File->matrix[irrep][0], coltot);

  return 0;

}
