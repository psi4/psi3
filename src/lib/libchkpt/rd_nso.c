/*!
  \file rd_nso.c
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

  psio_read_entry(PSIF_CHKPT, "::Num so", (char *) &nso, 
                  sizeof(int) );
  return nso;
}
