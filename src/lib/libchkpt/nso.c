/*!
  \file nso.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nso()  
** Reads in the total number of SOs.
**
** returns: int nso   total number of symmetry-adapted basis functions.
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
** arguments: 
**   \param int nso  total number of symmetry-adapted basis functions.
**
** returns: none
*/

void chkpt_wt_nso(int nso)
{
  psio_write_entry(PSIF_CHKPT, "::Num. SO", (char *) &nso, sizeof(int));
}
