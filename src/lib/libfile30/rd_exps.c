/*!
  \file rd_exps.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_exps():	Reads in the exponents of the primitive Gaussian functions.
**
**  takes no arguments.
**
**  returns: double *exps   The exponents are returned as an array of doubles.
*/


double *file30_rd_exps(void)
{
  double *exps;
  int nprim = 0;
  PSI_FPTR junk;
  PSI_FPTR exps_ptr;

  nprim = file30_rd_nprim();
  exps_ptr = (PSI_FPTR) (info30_.mpoint[4] - 1)*sizeof(int);

  exps = init_array(nprim);

  wreadw(info30_.filenum, (char *) exps, (int) nprim*sizeof(double),
	 exps_ptr, &junk);

  return exps;
}
