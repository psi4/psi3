#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file_mat_irrep_init(): Allocates and initializes memory for a
** matrix for a single irrep of a dpd two-electron file.
**
** Arguments:
**   struct dpdfile *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be prepared.
*/

int dpd_file_mat_irrep_init(struct dpdfile *File, int irrep)
{
  File->matrix[irrep] = block_matrix(File->params->rowtot[irrep],
				     File->params->coltot[irrep]);

  return 0;
}
