#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_buf_mat_irrep_init(): Allocates and initializes memory for a
** matrix for a single irrep of a dpd two-electron buffer.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
*/

int dpd_buf_mat_irrep_init(struct dpdbuf *Buf, int irrep)
{

  Buf->matrix[irrep] = block_matrix(Buf->params->rowtot[irrep],
				    Buf->params->coltot[irrep]);

  return 0;

}
