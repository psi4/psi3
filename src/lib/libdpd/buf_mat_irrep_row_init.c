#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_buf_mat_irrep_row_init(struct dpdbuf *Buf, int irrep)
{
  Buf->matrix[irrep] = block_matrix(1, Buf->params->coltot[irrep]);

  return 0;
}
