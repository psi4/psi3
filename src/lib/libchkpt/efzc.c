/*!
  \file efzc.c
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
*/
double chkpt_rd_efzc(void)
{
  double efzc;

  psio_read_entry(PSIF_CHKPT, "::Frozen core energy", (char *) &efzc, sizeof(double));

  return efzc;
}

/*!
** chkpt_wt_efzc(): Writes out the frozen-core energy.
**
** arguments: 
**  \param double efzc  the frozen-core energy.
**
** returns: none
*/
void chkpt_wt_efzc(double efzc)
{
  psio_write_entry(PSIF_CHKPT, "::Frozen core energy", (char *) &efzc, 
		   sizeof(double));
}

