#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "dpd.h"

int dpd_trans4_mat_irrep_shift31(dpdtrans4 *Trans, int irrep)
{
  int h, pq, Gr, Gs, r, nirreps;
  int rowtot, coltot;
  int *count;
  int *tempoff, *rowoff;
  double *data;

#ifdef DPD_TIMER
  timer_on("shift");
#endif
  if(Trans->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Trans->shift.shift_type);
      exit(Trans->shift.shift_type);
    }
  else Trans->shift.shift_type = 31;

  nirreps = Trans->buf.params->nirreps;
  rowtot = Trans->buf.params->coltot[irrep];
  coltot = Trans->buf.params->rowtot[irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Trans->matrix[irrep][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Trans->shift.coltot[irrep][h] = Trans->buf.params->qpi[h];
      Trans->shift.rowtot[irrep][h] = rowtot * Trans->buf.params->ppi[h^irrep];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Trans->shift.matrix[irrep] = (double ***) malloc(nirreps*sizeof(double **));
  for(h=0; h < nirreps; h++) 
      Trans->shift.matrix[irrep][h] =
	  ((!Trans->shift.rowtot[irrep][h]) ? NULL :
	   (double **) malloc(Trans->shift.rowtot[irrep][h] * sizeof(double *)));

  /* Calculate the row offsets */
  tempoff = init_int_array(nirreps);
  tempoff[0] = 0;
  for(h=1; h < nirreps; h++)
      tempoff[h] = tempoff[h-1] +
		   Trans->buf.params->ppi[h-1] *
		   Trans->buf.params->qpi[irrep^(h-1)];
  rowoff = init_int_array(nirreps);
  for(h=0; h < nirreps; h++)
      rowoff[h] = tempoff[h^irrep];
  free(tempoff);
  
  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over rows of original DPD matrix */
  for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {

      /* Loop over irreps of s */
      for(Gs=0; Gs < nirreps; Gs++) {
	  Gr = irrep^Gs;

	  /* Loop over orbitals in Gr */
	  for(r=0; (r < Trans->buf.params->ppi[Gr]) &&
		Trans->buf.params->qpi[Gs]; r++) {

	      /* Re-assign the row pointer */
              Trans->shift.matrix[irrep][Gs][count[Gs]] =
                &(data[pq*coltot + rowoff[Gs] +
		      (r * Trans->buf.params->qpi[Gs])]);

	      count[Gs]++;

	    }
	}
    }

  free(count); free(rowoff);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}
