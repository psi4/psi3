#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_file_mat_irrep_close_block(): Releases memory for a subblock of
** a matrix for a single irrep of a dpd two-electron file.
**
** Arguments:
**   struct dpdfile *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be freed.
**   int num_pq: The number of rows used.
**
*/


int dpd_file_mat_irrep_close_block(struct dpdfile *File, int irrep, int num_pq)
{
  if(num_pq && File->params->coltot[irrep]) free_block(File->matrix[irrep]);

  return 0;
}
