#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_buf_mat_irrep_init_block(): Allocates and initializes memory
** for a subblock of a matrix for a single irrep of a dpd two-electron
** buffer.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
**   int num_pq: The number of rows needed.
**
*/

int dpd_buf_mat_irrep_init_block(struct dpdbuf *Buf, int irrep, int num_pq)
{
  Buf->matrix[irrep] = block_matrix(num_pq,Buf->params->coltot[irrep]);

  return 0;

}
