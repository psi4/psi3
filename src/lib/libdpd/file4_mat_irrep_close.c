#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file4_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd four-index file.
**
** Arguments:
**   dpdfile4 *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be freed.
*/

int dpd_file4_mat_irrep_close(dpdfile4 *File, int irrep)
{
  int my_irrep;

  my_irrep = File->my_irrep;

  if(File->incore) return 0;  /* We need to keep the memory */

  if(File->params->rowtot[irrep] && File->params->coltot[irrep^my_irrep])
      free_block(File->matrix[irrep]);

  return 0;
}
