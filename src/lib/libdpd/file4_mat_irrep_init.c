#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file4_mat_irrep_init(): Allocates and initializes memory for a
** matrix for a single irrep of a dpd four-index file.
**
** Arguments:
**   dpdfile4 *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be prepared.
*/

int dpd_file4_mat_irrep_init(dpdfile4 *File, int irrep)
{
  int my_irrep;

  my_irrep = File->my_irrep;

  if(File->incore) return 0;  /* We've already got the memory */

  File->matrix[irrep] = block_matrix(File->params->rowtot[irrep],
				     File->params->coltot[irrep^my_irrep]);

  return 0;
}
