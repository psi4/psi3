#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

int dpd_buf4_mat_irrep_shift13(dpdbuf4 *Buf, int irrep)
{
  int h, i, nirreps;
  int *count, *dataoff;
  int rowtot, coltot;
  double *data;

#ifdef DPD_TIMER
  timer_on("shift");
#endif

  if(Buf->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Buf->shift.shift_type);
      exit(Buf->shift.shift_type);
    }
  else Buf->shift.shift_type = 13;

  nirreps = Buf->params->nirreps;
  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Buf->matrix[irrep][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Buf->shift.rowtot[irrep][h] = Buf->params->ppi[h];
      Buf->shift.coltot[irrep][h] = coltot * Buf->params->qpi[h^irrep];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Buf->shift.matrix[irrep] = (double ***) malloc(nirreps * sizeof(double **));
  for(h=0; h < nirreps; h++)
      Buf->shift.matrix[irrep][h] =
	   ((!Buf->shift.rowtot[irrep][h]) ? NULL :
	    (double **) malloc(Buf->shift.rowtot[irrep][h] * sizeof(double *)));

  /* Calculate the data offset */
  dataoff = init_int_array(nirreps);
  dataoff[0] = 0;
  for(h=1; h < nirreps; h++)
      dataoff[h] = dataoff[h-1] +
		   Buf->shift.rowtot[irrep][h-1] * Buf->shift.coltot[irrep][h-1];
		     

  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over irreps of isolated index */
  for(h=0; h < Buf->params->nirreps; h++) {
      for(i=0; (i < Buf->shift.rowtot[irrep][h]) && Buf->shift.coltot[irrep][h];
	  i++,count[h]++) {
          Buf->shift.matrix[irrep][h][count[h]] =
		    &(data[dataoff[h]+(Buf->shift.coltot[irrep][h])*i]);
        }
    }

  free(count); free(dataoff);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}
