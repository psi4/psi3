#include <stdio.h>
#include <string.h>
#include "dpd.h"

/* dpd_buf4_scmcopy(): Copies an existing four-index dpdbuf4 into another 
** file and multiplies it by a scalar at the same time.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the given dpd buffer.
**   int outfilenum: The PSI unit number for the new buffer.
**   char *label: A string labelling for this buffer.
**   double alpha: A scalar.
**
** NB: The buffer and file pq/rs parameters are assumed to be
** identical for the copy, obviously.  Hence, the anti flag must be off.
**
*/

int dpd_buf4_scmcopy(dpdbuf4 *InBuf, int outfilenum, char *label, double alpha)
{
  int h, row, col, my_irrep, rowtot, coltot;
  dpdbuf4 OutBuf;
  double *X;

  my_irrep = InBuf->file.my_irrep;

  dpd_buf4_init(&OutBuf, outfilenum, InBuf->file.my_irrep, InBuf->params->pqnum,
		InBuf->params->rsnum, InBuf->params->pqnum, 
                InBuf->params->rsnum, 0, label);

  for(h=0; h < InBuf->params->nirreps; h++) {

      dpd_buf4_mat_irrep_init(InBuf, h);
      dpd_buf4_mat_irrep_rd(InBuf, h);

      dpd_buf4_mat_irrep_init(&OutBuf, h);

      rowtot = InBuf->params->rowtot[h];
      coltot = InBuf->params->coltot[h^my_irrep];

      if(rowtot && coltot) {
          memcpy((void *) &(OutBuf.matrix[h][0][0]),
                 (const void *) &(InBuf->matrix[h][0][0]),
                 sizeof(double)*rowtot*coltot);
          C_DSCAL(rowtot*coltot, alpha, &(OutBuf.matrix[h][0][0]), 1);
       }

      dpd_buf4_mat_irrep_wrt(&OutBuf, h);

      dpd_buf4_mat_irrep_close(&OutBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }

  dpd_buf4_close(&OutBuf);

  return 0;
}
