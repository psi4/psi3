/*!
  \file rd_nmo.c
*/

#include "chkpt.h"
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nmo()  
** Reads in the total number of molecular orbitals.
**
** returns: int nmo   total number of molecular orbitals.
*/


int chkpt_rd_nmo(void)
{
  int nmo;

  psio_read_entry(PSIF_CHKPT, "::Num mo", (char *) &nmo, 
                  sizeof(int) );
  return nmo;
}
