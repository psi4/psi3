#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd two-electron file.
**
** Arguments:
**   struct dpdfile *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be freed.
*/

int dpd_file_mat_irrep_close(struct dpdfile *File, int irrep)
{
  if(File->params->rowtot[irrep] && File->params->coltot[irrep])
      free_block(File->matrix[irrep]);

  return 0;
}
