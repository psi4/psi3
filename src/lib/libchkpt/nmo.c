/*!
  \file nmo.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nmo()  
** Reads in the total number of molecular orbitals.
**
** returns: nmo = total number of molecular orbitals.
** \ingroup (CHKPT)
*/

int chkpt_rd_nmo(void)
{
  int nmo;
  char *keyword;
  keyword = chkpt_build_keyword("Num. MO's");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nmo, sizeof(int));

  free(keyword);
  return nmo;
}


/*!
** void chkpt_wt_nmo(int)  
** Writes out the total number of molecular orbitals.
**
** \param nmo = total number of molecular orbitals.
**
** \ingroup (CHKPT)
*/

void chkpt_wt_nmo(int nmo)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. MO's");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nmo, sizeof(int));

  free(keyword);
}
