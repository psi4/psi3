#include <stdio.h>
#include <libqt/qt.h>
#include "dpd.h"

int dpd_trans4_mat_irrep_wrt(dpdtrans4 *Trans, int irrep)
{
  int pq, rs;
  dpdbuf4 *Buf;

  Buf = &(Trans->buf);

#ifdef DPD_TIMER
  timer_on("trans4_rw");
#endif

  /* Loop over rows of transpose */
  for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
      for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
	  Buf->matrix[irrep][rs][pq] = Trans->matrix[irrep][pq][rs];
	}
    }

#ifdef DPD_TIMER
  timer_off("trans4_rw");
#endif

  return 0;
}
