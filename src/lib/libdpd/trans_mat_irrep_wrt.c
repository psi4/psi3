#include <stdio.h>
#include "dpd.h"

int dpd_trans_mat_irrep_wrt(struct dpdtrans *Trans, int irrep)
{
  int pq, rs;
  struct dpdbuf *Buf;

  Buf = &(Trans->buf);

  /* Loop over rows of transpose */
  for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
      for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
	  Buf->matrix[irrep][rs][pq] = Trans->matrix[irrep][pq][rs];
	}
    }

  return 0;
}
