#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_row_close(dpdbuf4 *Buf, int irrep)
{
  if(Buf->params->coltot[irrep])
    dpd_free_block(Buf->matrix[irrep], 1, Buf->params->coltot[irrep]);

  return 0;
}
