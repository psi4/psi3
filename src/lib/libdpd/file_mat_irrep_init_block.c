#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file_mat_irrep_init_block(): Allocates and initializes memory
** for a limited number of rows of amatrix for a single irrep of a dpd
** two-electron file.
**
** Arguments:
**   struct dpdfile *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be prepared.
**   int num_pq: The number of rows desired.
**
*/

int dpd_file_mat_irrep_init_block(struct dpdfile *File, int irrep, int num_pq)
{
  File->matrix[irrep] = block_matrix(num_pq,File->params->coltot[irrep]);

  return 0;
}
