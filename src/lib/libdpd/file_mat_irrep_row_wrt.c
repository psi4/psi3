#include <stdio.h>
#include <psio.h>
#include "dpd.h"

int dpd_file_mat_irrep_row_wrt(struct dpdfile *File, int irrep, int row)
{
  int coltot;
  psio_address irrep_ptr, row_ptr, next_address;
  
  irrep_ptr = File->lfiles[irrep];
  coltot = File->params->coltot[irrep];

  row_ptr = psio_get_address(irrep_ptr,row*coltot*sizeof(double));

  if(coltot) 
      psio_write(File->filenum, File->label, (char *) File->matrix[irrep][0],
		 coltot*sizeof(double), row_ptr, &next_address);

  return 0;

}
