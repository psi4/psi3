#include <stdio.h>
#include <qt.h>
#include "dpd.h"

int dpd_trans4_mat_irrep_wrt(dpdtrans4 *Trans, int irrep)
{
  int pq, rs;
  dpdbuf4 *Buf;

  Buf = &(Trans->buf);

  timer_on("trans4_rw");

  /* Loop over rows of transpose */
  for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
      for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
	  Buf->matrix[irrep][rs][pq] = Trans->matrix[irrep][pq][rs];
	}
    }

  timer_off("trans4_rw");

  return 0;
}
