#include <stdio.h>
#include "dpd.h"
#include <qt.h>

/* dpd_buf4_scm(): Multiplies every element of a four-index dpdbuf by a scalar.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the dpdbuf.
**   double alpha: The scalar.
**
** NB: This function is sometimes called automatically by the contractXXX
** functions to zero-out data that does not yet exist on disk.  In such a
** case, since each symmetry block is handled separately, it is possible 
** for only the first (totally symmetric) block to exist on disk while the
** others have yet to be created.  In such cases, the buf4_mat_irrep_rd()
** request will fail in libpsio, because the correct TOC entry exists
** (created when the first symmetry block was written), but the length of 
** the TOC entry will be too short for the next symmetry block to be added.
** So, to avoid this problem, here we manually check in the beginning to
** see if the TOC entry exists and should be read before the multiplication.
**
** TDC
** June 2000
**
*/

int dpd_buf4_scm(dpdbuf4 *InBuf, double alpha)
{
  int h, nirreps, my_irrep, new_buf4;
  int row, col, length;
  double *X;

  nirreps = InBuf->params->nirreps;
  my_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_scm");
#endif

  /* Look first for the TOC entry on disk */
  if(psio_tocscan(InBuf->file.filenum, InBuf->file.label) == NULL)
     new_buf4 = 1;
  else new_buf4 = 0;

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(InBuf, h);

      if(!new_buf4) dpd_buf4_mat_irrep_rd(InBuf, h);

      length = InBuf->params->rowtot[h] * InBuf->params->coltot[h];
      if(length) {
          X = &(InBuf->matrix[h][0][0]);
          C_DSCAL(length, alpha, X, 1);
       }

      dpd_buf4_mat_irrep_wrt(InBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }

#ifdef DPD_TIMER
  timer_off("buf4_scm");
#endif

  return 0;
}
