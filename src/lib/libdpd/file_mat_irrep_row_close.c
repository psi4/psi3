#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_file_mat_irrep_row_close(struct dpdfile *File, int irrep)
{
  if(File->params->coltot[irrep]) free_block(File->matrix[irrep]);

  return 0;
}
