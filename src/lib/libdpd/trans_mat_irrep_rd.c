#include <stdio.h>
#include "dpd.h"

int dpd_trans_mat_irrep_rd(struct dpdtrans *Trans, int irrep)
{
  int pq, rs;
  struct dpdbuf *Buf;

  Buf = &(Trans->buf);

  /* Loop over rows of transpose */
  for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
      for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
	  Trans->matrix[irrep][pq][rs] = Buf->matrix[irrep][rs][pq];
	}
    }

  return 0;
}
