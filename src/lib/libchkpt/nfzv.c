/*!
  \file nfzv.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nfzv()  
** Reads in the total number of frozen unoccupied molecular orbitals.
**
** returns: nfzv = total number of frozen unoccupied molecular orbitals.
** \ingroup (CHKPT)
*/

int chkpt_rd_nfzv(void)
{
  int nfzv;

  psio_read_entry(PSIF_CHKPT, "::Num. Frozen UOCC", (char *) &nfzv, 
                  sizeof(int) );
  return nfzv;
}


/*!
** void chkpt_wt_nfzv(int)  
** Writes out the total number of frozen unoccupied molecular orbitals.
**
** \param nfzv = total number of frozen unoccupied molecular orbitals.
**
** \ingroup (CHKPT)
*/

void chkpt_wt_nfzv(int nfzv)
{
  psio_write_entry(PSIF_CHKPT, "::Num. Frozen UOCC", (char *) &nfzv, sizeof(int));
}
