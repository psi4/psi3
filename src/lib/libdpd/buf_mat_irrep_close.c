#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_buf_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd two-electron buffer.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be freed.
**
** Note that shift information is freed here as well.
*/

int dpd_buf_mat_irrep_close(struct dpdbuf *Buf, int irrep)
{
  int h, nirreps;

  nirreps = Buf->params->nirreps;

  /* Free the shift structure for this irrep if used */
  if(Buf->shift.shift_type) {
      for(h=0; h < nirreps; h++) 
	  if(Buf->shift.rowtot[irrep][h])
	      free(Buf->shift.matrix[irrep][h]);
      free(Buf->shift.matrix[irrep]);
      Buf->shift.shift_type = 0;
    }

  if(Buf->params->rowtot[irrep] && Buf->params->coltot[irrep])
      free_block(Buf->matrix[irrep]);

  return 0;
}
