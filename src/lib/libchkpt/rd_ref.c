/*!
  \file rd_ref.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_ref()  
** Reads the reference type from the flag in checkpoint
** 0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF 
**
** returns: int refnum   number indicating the reference.
*/

int chkpt_rd_ref(void)
{
  int refnum;

  psio_read_entry(PSIF_CHKPT, "::Reference", (char *) &refnum, 
                  sizeof(int) );
  return refnum;
}
