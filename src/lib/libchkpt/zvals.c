/*!
    \file zvals.c
*/

#include <stdio.h>
#include <stdlib.h>
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
  char *keyword;
  keyword = chkpt_build_keyword("Nuclear charges");

  natom = chkpt_rd_natom();
  zvals = init_array(natom);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) zvals, 
    natom*sizeof(double));

  free(keyword);
  return zvals;
}

/*!
** chkpt_wt_zvals()
** Writes the nuclear charges to the checkpoint file.
**
** \param double *zvals: An array of the charges
**
** returns: nothing
*/

void chkpt_wt_zvals(double *zvals)
{
  int natom;
  char *keyword;
  keyword = chkpt_build_keyword("Nuclear charges");

  natom = chkpt_rd_natom();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) zvals, 
    natom*sizeof(double));

  free(keyword);
}
