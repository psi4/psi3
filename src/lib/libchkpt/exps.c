/*!
  \file exps.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>

/*!
** chkpt_rd_exps():	
** Reads in the exponents of the primitive Gaussian functions.
**
** takes no arguments.
**
** returns: double *exps   
** The exponents are returned as an array of doubles.
** \ingroup(CHKPT)
*/

double *chkpt_rd_exps(void)
{
  double *exps;
  int nprim = 0;

  nprim = chkpt_rd_nprim();

  exps = init_array(nprim);

  psio_read_entry(PSIF_CHKPT, "::Exponents", (char *) exps, 
		  nprim*sizeof(double));

  return exps;
}


/*!
** chkpt_wt_exps(): 
** Writes out the exponents of the primitive Gaussian functions.
**
** arguments:
**  \param exps = The exponents are returned as an array of doubles.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_exps(double *exps)
{
  int nprim;

  nprim = chkpt_rd_nprim();

  psio_write_entry(PSIF_CHKPT, "::Exponents", (char *) exps, 
		   nprim*sizeof(double));
}
