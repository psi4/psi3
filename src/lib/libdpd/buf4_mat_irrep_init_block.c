#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_buf4_mat_irrep_init_block(): Allocates and initializes memory
** for a subblock of a matrix for a single irrep of a dpd four-index
** buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
**   int num_pq: The number of rows needed.
**
*/

int dpd_buf4_mat_irrep_init_block(dpdbuf4 *Buf, int irrep, int num_pq)
{
  Buf->matrix[irrep] = dpd_block_matrix(num_pq,Buf->params->coltot[irrep]);

  return 0;

}
