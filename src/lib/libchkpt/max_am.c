/*!
  \file max_am.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_max_am()
** Reads in the maximum orbital quantum number 
** of AOs in the basis.
**
** returns: int max_am (0 corresponds to s-functions, 
**                      1 - to up to p-functions, etc.)
*/


int chkpt_rd_max_am(void)
{
  int max_am;

  psio_read_entry(PSIF_CHKPT, "::Max. AM", (char *) &max_am, sizeof(int));
  return max_am;
}

/*!
** void chkpt_wt_max_am()
** Writes out the maximum orbital quantum number 
** of AOs in the basis.
**
** arguments: 
**  \param int max_am (0 corresponds to s-functions, 
**                     1 - to up to p-functions, etc.)
**
** returns: none
*/


void chkpt_wt_max_am(int max_am)
{
  psio_write_entry(PSIF_CHKPT, "::Max. AM", (char *) &max_am, sizeof(int));
}
