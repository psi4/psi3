/*!
  \file calcs.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_ncalcs()  
** Reads in the total number of calculations.
**
** returns: ncalcs = total number of calculations in checkpoint
** \ingroup (CHKPT)
*/

int chkpt_rd_ncalcs(void)
{
  int ncalcs;

  psio_read_entry(PSIF_CHKPT, "::Num calcs", (char *) &ncalcs, 
                  sizeof(int) );
  return ncalcs;
}
