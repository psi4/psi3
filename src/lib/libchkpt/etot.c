/*!
  \file etot.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_etot(): Reads in the total energy.
**
**  takes no arguments.
**
**  returns: double etot  the total energy.
**  \ingroup (CHKPT)
*/

double chkpt_rd_etot(void)
{
  double etot;
  char *keyword;
  keyword = chkpt_build_keyword("Total energy");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

  free(keyword);
  return etot;
}

/*!
** chkpt_wt_etot(): Writes out the total energy.
**
**  arguments: 
**   \param double etot  the total energy.
**
**  returns: none
**  \ingroup (CHKPT)
*/

void chkpt_wt_etot(double etot)
{
  char *keyword;
  keyword = chkpt_build_keyword("Total energy");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

  free(keyword);
}
