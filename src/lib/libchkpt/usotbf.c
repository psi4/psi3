/*!
  \file usotbf.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_usotbf(): Reads in the SO to basis functions transformation matrix 
**
** takes no arguments.
**
** returns: usotbf = Read in a num_so by num_so matrix of doubles
**
** \ingroup (CHKPT)
*/

double **chkpt_rd_usotbf(void)
{
  double **usotbf;
  int num_so, i;
  psio_address ptr;

  num_so = chkpt_rd_nso();

  usotbf = block_matrix(num_so,num_so);

  ptr = PSIO_ZERO;
  for(i=0;i<num_so;i++)
    psio_read(PSIF_CHKPT, "::SO->BF transmat", (char *) usotbf[i], 
	      (int) num_so*sizeof(double), ptr, &ptr);

  return usotbf;
}


/*!
** chkpt_wt_usotbf(): Writes out the SO to basis functions transformation 
**                    matrix 
**
** \param usotbf = A num_so by num_so matrix of doubles
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_usotbf(double **usotbf)
{
  int num_so, i;
  psio_address ptr;

  num_so = chkpt_rd_nso();

  ptr = PSIO_ZERO;
  for(i=0;i<num_so;i++)
    psio_write(PSIF_CHKPT, "::SO->BF transmat", (char *) usotbf[i], 
	       (int) num_so*sizeof(double), ptr, &ptr);
}
