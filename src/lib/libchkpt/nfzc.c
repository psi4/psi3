/*!
  \file nfzc.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nfzc()  
** Reads in the total number of frozen doubly occupied molecular orbitals.
**
** returns: nfzc = total number of frozen doubly occupied molecular orbitals.
** \ingroup (CHKPT)
*/

int chkpt_rd_nfzc(void)
{
  int nfzc;

  psio_read_entry(PSIF_CHKPT, "::Num. Frozen DOCC", (char *) &nfzc, 
                  sizeof(int) );
  return nfzc;
}


/*!
** void chkpt_wt_nfzc(int)  
** Writes out the total number of frozen doubly occupied molecular orbitals.
**
** \param nfzc = total number of frozen doubly occupied molecular orbitals.
**
** \ingroup (CHKPT)
*/

void chkpt_wt_nfzc(int nfzc)
{
  psio_write_entry(PSIF_CHKPT, "::Num. Frozen DOCC", (char *) &nfzc, sizeof(int));
}
