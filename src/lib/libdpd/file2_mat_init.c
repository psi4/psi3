#include "dpd.h"

int dpd_file2_mat_init(dpdfile2 *File)
{
  int h, my_irrep;

  my_irrep = File->my_irrep;

  if(File->incore) return 0;  /* We've already got the memory */

  for(h=0; h < File->params->nirreps; h++)
      File->matrix[h] = block_matrix(File->params->rowtot[h],
				     File->params->coltot[h^my_irrep]);

  return 0;
}
