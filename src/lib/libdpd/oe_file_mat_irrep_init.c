#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_oe_file_mat_irrep_init(struct oe_dpdfile *File, int irrep)
{
  File->matrix[irrep] = block_matrix(File->params->rowtot[irrep],
				     File->params->coltot[irrep]);

  return 0;
}
