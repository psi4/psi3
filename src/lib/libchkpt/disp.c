/*!
  \file disp.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_disp():  Reads in the current geometry displacement number.
**
**   takes no arguments.
**
**   returns: int disp   the current geometry displacement number
*/
int chkpt_rd_disp(void)
{
  int disp;

  psio_read_entry(PSIF_CHKPT, "::Current displacement", (char *) &disp,
		  sizeof(int));
  return disp;
}

/*!
** chkpt_wt_disp():  Writes out the current geometry displacement number.
**
**  arguments: 
**   \param int disp   the current geometry displacement number
**
** returns: none
*/
void chkpt_wt_disp(int disp)
{
  psio_write_entry(PSIF_CHKPT, "::Current displacement", (char *) &disp,
		   sizeof(int));
}
