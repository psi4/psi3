#include "dpd.h"

int dpd_file2_mat_wrt(dpdfile2 *File)
{
  int h, my_irrep, rowtot, coltot;
  psio_address irrep_ptr, next_address;

  my_irrep = File->my_irrep;

  if(File->incore) return 0;  /* We're keeping this data in core */

  for(h=0; h < File->params->nirreps; h++) {
      irrep_ptr = File->lfiles[h];
      rowtot = File->params->rowtot[h];
      coltot = File->params->coltot[h^my_irrep];

      if(rowtot && coltot)
	  psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
		     rowtot*coltot*sizeof(double), irrep_ptr, &next_address);
    }

  return 0;
}
