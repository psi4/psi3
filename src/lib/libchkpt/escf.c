/*!
  \file escf.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_escf(): Reads in the scf energy.
**
**  takes no arguments.
**
**  returns: double escf  the scf energy.
** \ingroup (CHKPT)
*/

double chkpt_rd_escf(void)
{
  double escf;
  char *keyword;
  keyword = chkpt_build_keyword("SCF energy");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

  free(keyword);
  return escf;
}

/*!
** chkpt_wt_escf(): Writes out the scf energy.
**
**  arguments: 
**   \param double escf  the scf energy.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_escf(double escf)
{
  char *keyword;
  keyword = chkpt_build_keyword("SCF energy");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

  free(keyword);
}
