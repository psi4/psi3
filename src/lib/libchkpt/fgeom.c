/*!
  \file fgeom.c
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_fgeom():  Reads in full cartesian geometry including dummy atoms
**
**   takes no arguments.
**   returns: double **full_geom;
**     
*/

double **chkpt_rd_fgeom(void)
{
  int nentry;
  double **fgeom;

  nentry = chkpt_rd_nentry();

  fgeom = block_matrix(nentry,3);

  psio_read_entry(PSIF_CHKPT, "::Full cartesian geometry", (char *) fgeom[0], 
		  (int) 3*nentry*sizeof(double));

  return  fgeom;
}

/*!
** chkpt_wt_fgeom():  Writes out full cartesian geometry including dummy atoms
**
**  arguments: 
**   \param double **full_geom;
**
** returns: none
**     
*/

void chkpt_wt_fgeom(double **fgeom)
{
  int nentry;

  nentry = chkpt_rd_nentry();

  psio_write_entry(PSIF_CHKPT, "::Full cartesian geometry", (char *) fgeom[0], 
		  (int) 3*nentry*sizeof(double));
}
