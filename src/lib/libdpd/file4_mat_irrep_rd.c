#include <stdio.h>
#include <psio.h>
#include <qt.h>
#include "dpd.h"

int dpd_file4_mat_irrep_rd(dpdfile4 *File, int irrep)
{
  int rowtot, coltot, my_irrep;
  psio_address irrep_ptr, next_address;

  if(File->incore) return 0;  /* We already have this data in core */

  /* If the data doesn't actually exist on disk, we just leave */
  if(psio_tocscan(File->filenum, File->label) == NULL) return 1;

  timer_on("file4_rd");

  my_irrep = File->my_irrep;
  irrep_ptr = File->lfiles[irrep];
  rowtot = File->params->rowtot[irrep];
  coltot = File->params->coltot[irrep^my_irrep];

  if(rowtot && coltot)
     psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       rowtot*coltot*sizeof(double), irrep_ptr, &next_address);

  timer_off("file4_rd");

  return 0;

}
