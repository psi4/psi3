#include <stdio.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

/* dpd_buf4_mat_irrep_init(): Allocates and initializes memory for a
** matrix for a single irrep of a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
*/

int dpd_buf4_mat_irrep_init(dpdbuf4 *Buf, int irrep)
{
  int my_irrep, rowtot, coltot;

  my_irrep = Buf->file.my_irrep;
  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep^my_irrep];

  timer_on("buf4_init");

  if(rowtot*coltot) {

      /* If the file member is already in cache and its ordering is the 
         same as the parent buffer, don't malloc() memory, just assign 
         the pointer */
      if(Buf->file.incore && !(Buf->anti) && 
          (Buf->params->pqnum == Buf->file.params->pqnum) &&
          (Buf->params->rsnum == Buf->file.params->rsnum))
          Buf->matrix[irrep] = Buf->file.matrix[irrep];
      else 
          Buf->matrix[irrep] = block_matrix(rowtot,coltot);

    }

  timer_off("buf4_init");

  return 0;

}
