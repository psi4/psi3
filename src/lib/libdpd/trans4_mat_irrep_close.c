#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_trans4_mat_irrep_close(dpdtrans4 *Trans, int irrep)
{
  int h, nirreps, rowtot, coltot;

  nirreps = Trans->buf.params->nirreps;
  rowtot = Trans->buf.params->coltot[irrep];
  coltot = Trans->buf.params->rowtot[irrep];

  /* Free the shift structure for this irrep if used */
  if(Trans->shift.shift_type) {
      for(h=0; h < nirreps; h++)
	  if(Trans->shift.rowtot[irrep][h])
	      free(Trans->shift.matrix[irrep][h]);
      free(Trans->shift.matrix[irrep]);
      Trans->shift.shift_type = 0;
    }

  if(rowtot * coltot)
      dpd_free_block(Trans->matrix[irrep], rowtot, coltot);
  
  return 0;
}
