/*!
  \file fgeom.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_fgeom():  Reads in full cartesian geometry including dummy atoms
**
** takes no arguments.
** returns: double **full_geom;
** \ingroup (CHKPT)    
*/

double **chkpt_rd_fgeom(void)
{
  int nallatom;
  double **fgeom;
  char keyword[PSIO_KEYLEN];

  nallatom = chkpt_rd_nallatom();

  fgeom = block_matrix(nallatom,3);

  /* build keyword string */
  sprintf(keyword, ":%s:%s", chkpt_prefix, "Full cartesian geometry");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) fgeom[0], (int) 3*nallatom*sizeof(double));

  return  fgeom;
}

/*!
** chkpt_wt_fgeom():  Writes out full cartesian geometry including dummy atoms
**
** arguments: 
**   \param full_geom = Matrix for cartesian coordinates
**
** returns: none
** \ingroup (CHKPT) 
*/

void chkpt_wt_fgeom(double **fgeom)
{
  int nallatom;
  char keyword[PSIO_KEYLEN];

  nallatom = chkpt_rd_nallatom();

  sprintf(keyword, ":%s:%s", chkpt_prefix, "Full cartesian geometry");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) fgeom[0], (int) 3*nallatom*sizeof(double));
}
