#include <stdio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_row_init(dpdbuf4 *Buf, int irrep)
{
#ifdef DPD_TIMER
  timer_on("b4_rowinit");
#endif
  Buf->matrix[irrep] = dpd_block_matrix(1, Buf->params->coltot[irrep]);
#ifdef DPD_TIMER
  timer_off("b4_rowinit");
#endif

  return 0;
}
