#include <stdio.h>
#include <psio.h>
#include <qt.h>
#include "dpd.h"

int dpd_file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row)
{
  int coltot, my_irrep;
  psio_address irrep_ptr, row_ptr, next_address;

  if(File->incore) return 0;  /* We already have this data in core */

#ifdef DPD_TIMER
  timer_on("f4_rowrd");
#endif

  my_irrep = File->my_irrep;

  irrep_ptr = File->lfiles[irrep];
  coltot = File->params->coltot[irrep^my_irrep];

  row_ptr = psio_get_address(irrep_ptr, row*coltot*sizeof(double));

  if(coltot) 
      psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
		coltot*sizeof(double), row_ptr, &next_address);

#ifdef DPD_TIMER
  timer_off("f4_rowrd");
#endif

  return 0;

}
