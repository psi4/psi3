/*!
  \file wt_zvals.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_zvals():  Writes in an array of nuclear charges.
**
**  arguments 
** \param double *zvals	An array of natom doubles which contains
**			nuclear charges.
**
**   returns no arguments.
*/


void file30_wt_zvals(double *zvals)
{
  int natom;
  PSI_FPTR zvals_ptr, junk;

  zvals_ptr = (PSI_FPTR) (info30_.mpoint[0] - 1)*sizeof(int);

  wwritw(info30_.filenum, (char *) zvals, (int) sizeof(double)*file30_rd_natom(),
         zvals_ptr, &junk);
}
