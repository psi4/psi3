#include <stdio.h>
#include <string.h>
#include "dpd.h"

/* dpd_buf4_copy(): Copies an existing four-index dpdbuf4 into another file.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the given dpd buffer.
**   int outfilenum: The PSI unit number for the new buffer.
**   char *label: A string labelling for this buffer.
**
** NB: The buffer and file pq/rs parameters are assumed to be
** identical for the copy, obviously.  Hence, anti flag must be off.
**
** Converted to use buf4 only rather than assumptions about file4.
** TDC, September 1999
*/

int dpd_buf4_copy(dpdbuf4 *InBuf, int outfilenum, char *label)
{
  int h, row, col, my_irrep;
  long int rowtot, coltot;
  dpdbuf4 OutBuf;

  my_irrep = InBuf->file.my_irrep;

  dpd_buf4_init(&OutBuf, outfilenum, InBuf->file.my_irrep, InBuf->params->pqnum,
		InBuf->params->rsnum, InBuf->params->pqnum, 
                InBuf->params->rsnum, 0, label);

  for(h=0; h < InBuf->params->nirreps; h++) {

      dpd_buf4_mat_irrep_init(InBuf, h);
      dpd_buf4_mat_irrep_rd(InBuf, h);

      dpd_buf4_mat_irrep_init(&OutBuf, h);

/*
      for(row=0; row < InBuf->params->rowtot[h]; row++) {
	  for(col=0; col < InBuf->params->coltot[h^my_irrep]; col++) {
	      OutBuf.matrix[h][row][col] = InBuf->matrix[h][row][col];
	    }
	}
*/

      /* Use memcpy() instead of element by element copy */
      rowtot = InBuf->params->rowtot[h];
      coltot = InBuf->params->coltot[h^my_irrep];

      if(rowtot && coltot) 
          memcpy((void *) &(OutBuf.matrix[h][0][0]),
                 (const void *) &(InBuf->matrix[h][0][0]),
                 sizeof(double)*rowtot*coltot);

      dpd_buf4_mat_irrep_wrt(&OutBuf, h);

      dpd_buf4_mat_irrep_close(&OutBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }

  dpd_buf4_close(&OutBuf);

  return 0;
}
