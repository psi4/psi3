#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <qt.h>
#include "dpd.h"

int dpd_trans4_mat_irrep_shift13(dpdtrans4 *Trans, int irrep)
{
  int h, i, nirreps;
  int *count, *dataoff;
  int rowtot, coltot;
  double *data;

#ifdef DPD_TIMER
  timer_on("shift");
#endif

  if(Trans->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Trans->shift.shift_type);
      exit(Trans->shift.shift_type);
    }
  else Trans->shift.shift_type = 13;

  nirreps = Trans->buf.params->nirreps;
  rowtot = Trans->buf.params->coltot[irrep];
  coltot = Trans->buf.params->rowtot[irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Trans->matrix[irrep][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Trans->shift.rowtot[irrep][h] = Trans->buf.params->rpi[h];
      Trans->shift.coltot[irrep][h] = coltot * Trans->buf.params->spi[h^irrep];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Trans->shift.matrix[irrep] = (double ***) malloc(nirreps * sizeof(double **));
  for(h=0; h < nirreps; h++)
      Trans->shift.matrix[irrep][h] =
	   ((!Trans->shift.rowtot[irrep][h]) ? NULL :
	    (double **) malloc(Trans->shift.rowtot[irrep][h] * sizeof(double *)));

  /* Calculate the data offset */
  dataoff = init_int_array(nirreps);
  dataoff[0] = 0;
  for(h=1; h < nirreps; h++)
      dataoff[h] = dataoff[h-1] +
		   Trans->shift.rowtot[irrep][h-1] *
		   Trans->shift.coltot[irrep][h-1];
		     

  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over irreps of isolated index */
  for(h=0; h < nirreps; h++) {
      for(i=0; (i < Trans->shift.rowtot[irrep][h]) &&
	    Trans->shift.coltot[irrep][h]; i++,count[h]++) {
          Trans->shift.matrix[irrep][h][count[h]] =
		    &(data[dataoff[h]+(Trans->shift.coltot[irrep][h])*i]);
        }
    }

  free(count); free(dataoff);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}
