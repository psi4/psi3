#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

int dpd_file4_mat_irrep_wrt_block(dpdfile4 *File, int irrep, int start_pq,
				 int num_pq)
{
  int rowtot, coltot, my_irrep;
  psio_address irrep_ptr, next_address;

  if(File->incore) {
      dpd_file4_cache_dirty(File);  /* Flag this cache entry for writing */
      return 0;  /* We're keeping this data in core */
    }

  my_irrep = File->my_irrep;
  irrep_ptr = File->lfiles[irrep];
  rowtot = num_pq;
  coltot = File->params->coltot[irrep^my_irrep];

  /* Advance file pointer to current row */
  irrep_ptr = psio_get_address(irrep_ptr, start_pq*coltot*sizeof(double));

  if(rowtot && coltot)
     psio_write(File->filenum, File->label, (char *) File->matrix[irrep][0],
		rowtot*coltot*sizeof(double), irrep_ptr, &next_address);

  return 0;

}
