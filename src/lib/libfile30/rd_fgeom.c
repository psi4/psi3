/*!
  \file rd_fgeom.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_fgeom():  Reads in full cartesian geometry including dummy atoms
**
**   takes no arguments.
**   returns: double **full_geom;
**     
*/

 double **file30_rd_fgeom(void)
{
  int nentry;
  PSI_FPTR fgeom_ptr, junk;
  double **fgeom;

  nentry = file30_rd_nentry();

  fgeom = block_matrix(nentry,3);

  fgeom_ptr = (PSI_FPTR) (info30_.mpoint[50] - 1)*sizeof(int);

  wreadw(info30_.filenum, (char *) fgeom[0], (int) 3*nentry*sizeof(double),
           fgeom_ptr, &junk);

  return  fgeom;
}
