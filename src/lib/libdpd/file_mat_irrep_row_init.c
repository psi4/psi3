#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_file_mat_irrep_row_init(struct dpdfile *File, int irrep)
{
  File->matrix[irrep] = block_matrix(1, File->params->coltot[irrep]);

  return 0;
}
