/*!
  \file rd_phase_check.c
*/

#include "chkpt.h"
#include <libpsio/psio.h>

/*!
** int chkpt_rd_phase_check()
**
** THIS FUNCTION NEEDS DOCUMENTATION!!!
** VMR 23/Mar/2002
*/

int chkpt_rd_phase_check(void)
{
  int pcheck;

  psio_read_entry(PSIF_CHKPT, "::Phase check", (char *) &pcheck, 
                  sizeof(int) );
  return pcheck;
}
