#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_row_zero(dpdbuf4 *Buf, int irrep, int row)
{
  int coltot;
  
  coltot = Buf->params->coltot[irrep];

  if(coltot)
      zero_arr(Buf->matrix[irrep][0], coltot);

  return 0;

}
