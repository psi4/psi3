/*!
  \file nso.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nso()  
** Reads in the total number of SOs.
**
** returns: nso = total number of symmetry-adapted basis functions.
** \ingroup (CHKPT)
*/

int chkpt_rd_nso(void)
{
  int nso;
  char *keyword;
  keyword = chkpt_build_keyword("Num. SO");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nso, sizeof(int));

  free(keyword);
  return nso;
}


/*!
** void chkpt_wt_nso(int)  
** Writes out the total number of SOs.
**
** \param nso = total number of symmetry-adapted basis functions.
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_nso(int nso)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. SO");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nso, sizeof(int));

  free(keyword);
}
