/*!
  \file rd_max_am.c
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

  psio_read_entry(PSIF_CHKPT, "::Max am", (char *) &max_am, 
                  sizeof(int) );
  return max_am;
}
