#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_oe_file_mat_irrep_close(struct oe_dpdfile *File, int irrep)
{
  if(File->params->rowtot[irrep] && File->params->coltot[irrep])
      free_block(File->matrix[irrep]);

  return 0;
}
