#include <stdio.h>
#include <psio.h>
#include "dpd.h"

int dpd_file4_mat_irrep_rd_block(dpdfile4 *File, int irrep, int start_pq,
				int num_pq)
{
  int rowtot, coltot;
  psio_address irrep_ptr, next_address;

  if(File->incore) return 0;  /* We already have this data in core */

  irrep_ptr = File->lfiles[irrep];
  rowtot = num_pq;
  coltot = File->params->coltot[irrep];

  /* Advance file pointer to current row */
  irrep_ptr = psio_get_address(irrep_ptr, start_pq*coltot*sizeof(double));

  if(rowtot && coltot)
     psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       rowtot*coltot*sizeof(double), irrep_ptr, &next_address);

  return 0;

}
