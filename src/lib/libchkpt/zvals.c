/*!
    \file zvals.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_zvals()
** Reads the nuclear charges from the checkpoint file.
**
** arguments: none
**
** returns: 
**   double *zvals: An array of the charges
*/

double *chkpt_rd_zvals(void)
{
  int natom;
  double *zvals;

  natom = chkpt_rd_natom();

  zvals = init_array(natom);

  psio_read_entry(PSIF_CHKPT, "::Nuclear charges", (char *) zvals, natom*sizeof(double));

  return zvals;
}

/*!
** chkpt_wt_zvals()
** Writes the nuclear charges to the checkpoint file.
**
** arguments:
** \param double *zvals: An array of the charges
**
** returns: nothing
*/

void chkpt_wt_zvals(double *zvals)
{
  int natom;

  natom = chkpt_rd_natom();

  psio_write_entry(PSIF_CHKPT, "::Nuclear charges", (char *) zvals, natom*sizeof(double));
}
