/*!
  \file efzc.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_efzc(): Reads in the frozen-core energy.
**
**   takes no arguments.
**
**   returns: double efzc  the frozen-core energy.
** \ingroup (CHKPT)
*/
double chkpt_rd_efzc(void)
{
  double efzc;
  char *keyword;
  keyword = chkpt_build_keyword("Frozen core energy");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &efzc, sizeof(double));

  free(keyword);
  return efzc;
}


/*!
** chkpt_wt_efzc(): Writes out the frozen-core energy.
**
** \param efzc = the frozen-core energy.
**
** returns: none
** \ingroup (CHKPT)
*/
void chkpt_wt_efzc(double efzc)
{
  char *keyword;
  keyword = chkpt_build_keyword("Frozen core energy");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &efzc, 
    sizeof(double));

  free(keyword);
}

