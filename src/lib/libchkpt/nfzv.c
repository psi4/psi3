/*!
  \file nfzv.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. Frozen UOCC");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int) );

  free(keyword);
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. Frozen UOCC");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int));

  free(keyword);
}
