/*!
  \file efzc.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
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

  psio_read_entry(PSIF_CHKPT, "::Frozen core energy", (char *) &efzc, 
                  sizeof(double));

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
  psio_write_entry(PSIF_CHKPT, "::Frozen core energy", (char *) &efzc, 
		   sizeof(double));
}

