#include <stdio.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

int dpd_file4_mat_irrep_row_init(dpdfile4 *File, int irrep)
{
  int my_irrep;

  if(File->incore) return 0;  /* We already have the whole matrix in core */

  timer_on("f4_rowinit");

  my_irrep = File->my_irrep;
  
  File->matrix[irrep] = block_matrix(1, File->params->coltot[irrep^my_irrep]);

  timer_off("f4_rowinit");

  return 0;
}
