#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_buf_mat_irrep_row_close(struct dpdbuf *Buf, int irrep)
{
  if(Buf->params->coltot[irrep]) free_block(Buf->matrix[irrep]);

  return 0;
}
