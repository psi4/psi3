/*!
  \file nmo.c
*/

#include "chkpt.h"
#include <psifiles.h>
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

  psio_read_entry(PSIF_CHKPT, "::Num. MO's", (char *) &nmo, 
                  sizeof(int) );
  return nmo;
}

/*!
** void chkpt_wt_nmo(int)  
** Writes out the total number of molecular orbitals.
**
** arguments: 
**  \param int nmo   total number of molecular orbitals.
*/

void chkpt_wt_nmo(int nmo)
{
  psio_write_entry(PSIF_CHKPT, "::Num. MO's", (char *) &nmo, sizeof(int));
}
