#include <stdio.h>
#include <psio.h>
#include "dpd.h"

int dpd_file4_mat_irrep_row_wrt(dpdfile4 *File, int irrep, int row)
{
  int coltot, my_irrep;
  psio_address irrep_ptr, row_ptr, next_address;

  if(File->incore) {
      dpd_file4_cache_dirty(File);  /* Flag this cache entry for writing */
      return 0;  /* We're keeping the data in core */
    }

  my_irrep = File->my_irrep;
  
  irrep_ptr = File->lfiles[irrep];
  coltot = File->params->coltot[irrep^my_irrep];

  row_ptr = psio_get_address(irrep_ptr,row*coltot*sizeof(double));

  if(coltot) 
      psio_write(File->filenum, File->label, (char *) File->matrix[irrep][0],
		 coltot*sizeof(double), row_ptr, &next_address);

  return 0;

}
