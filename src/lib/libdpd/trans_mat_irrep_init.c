#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_trans_mat_irrep_init(struct dpdtrans *Trans, int irrep)
{
  Trans->matrix[irrep] = block_matrix(Trans->buf.params->coltot[irrep],
				      Trans->buf.params->rowtot[irrep]);

  return 0;
}
