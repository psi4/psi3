/*!
  \file rd_zvals.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_zvals():  Reads in the nuclear charge for each atom.
**
**   takes no arguments.
**
**   returns: double *zvals   An array natom long which contains the
**     nuclear charge (as a double) for each atom.
*/


double *file30_rd_zvals(void)
{
  int natom;
  PSI_FPTR zvals_ptr, junk;
  double *zvals;

  natom = file30_rd_natom();

  zvals_ptr = (PSI_FPTR) (info30_.mpoint[0] - 1)*sizeof(int);

  zvals = init_array(natom);

  wreadw(info30_.filenum, (char *) zvals, (int) sizeof(double)*natom,
         zvals_ptr, &junk);

  return zvals;
}
