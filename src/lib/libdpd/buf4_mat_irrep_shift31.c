#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_shift31(dpdbuf4 *Buf, int irrep)
{
  int h, pq, Gr, Gs, r, nirreps;
  int rowtot, coltot;
  int *count;
  int *tempoff, *rowoff;
  double *data;

  timer_on("shift");

  if(Buf->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Buf->shift.shift_type);
      exit(Buf->shift.shift_type);
    }
  else Buf->shift.shift_type = 31;

  nirreps = Buf->params->nirreps;
  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Buf->matrix[irrep][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Buf->shift.coltot[irrep][h] = Buf->params->spi[h];
      Buf->shift.rowtot[irrep][h] = rowtot * Buf->params->rpi[h^irrep];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Buf->shift.matrix[irrep] = (double ***) malloc(nirreps*sizeof(double **));
  for(h=0; h < nirreps; h++) 
      Buf->shift.matrix[irrep][h] =
	   ((!Buf->shift.rowtot[irrep][h]) ? NULL :
	    (double **) malloc(Buf->shift.rowtot[irrep][h] * sizeof(double *)));

  /* Calculate the row offsets */
  tempoff = init_int_array(nirreps);
  tempoff[0] = 0;
  for(h=1; h < nirreps; h++)
      tempoff[h] = tempoff[h-1] +
		     Buf->params->rpi[h-1] * Buf->params->spi[irrep^(h-1)];
  rowoff = init_int_array(nirreps);
  for(h=0; h < nirreps; h++)
      rowoff[h] = tempoff[h^irrep];
  free(tempoff);
  
  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over rows of original DPD matrix */
  for(pq=0; pq < Buf->params->rowtot[irrep]; pq++) {

      /* Loop over irreps of s */
      for(Gs=0; Gs < nirreps; Gs++) {
	  Gr = irrep^Gs;

	  /* Loop over orbitals in Gr */
	  for(r=0; (r < Buf->params->rpi[Gr]) && Buf->params->spi[Gs]; r++) {

	      /* Re-assign the row pointer */
              Buf->shift.matrix[irrep][Gs][count[Gs]] =
                &(data[pq*coltot + rowoff[Gs] + (r * Buf->params->spi[Gs])]);

	      count[Gs]++;

	    }
	}
    }

  free(count); free(rowoff);

  timer_off("shift");

  return 0;
}
