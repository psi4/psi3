#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_row_close(dpdbuf4 *Buf, int irrep)
{
  if(Buf->params->coltot[irrep]) free_block(Buf->matrix[irrep]);

  return 0;
}
