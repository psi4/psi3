/*!
  \file eref.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_eref(): Reads in the reference energy.
**
**   takes no arguments.
**
**   returns: double eref  the reference energy.
**
** \ingroup (CHKPT)
*/
double chkpt_rd_eref(void)
{
  double eref;
  char *keyword;
  keyword = chkpt_build_keyword("Reference energy");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &eref, sizeof(double));

  free(keyword);
  return eref;
}


/*!
** chkpt_wt_eref(): Writes out the reference energy.
**
** \param double eref = the reference energy.
**
** returns: none
**
** \ingroup (CHKPT)
*/
void chkpt_wt_eref(double eref)
{
  char *keyword;
  keyword = chkpt_build_keyword("Reference energy");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &eref, 
                   sizeof(double));

  free(keyword);
}
