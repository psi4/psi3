/*!
  \file nso.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
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

  psio_read_entry(PSIF_CHKPT, "::Num. SO", (char *) &nso, sizeof(int));
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
  psio_write_entry(PSIF_CHKPT, "::Num. SO", (char *) &nso, sizeof(int));
}
