/*!
  \file ref.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_ref()  
** Reads the reference type from the flag in checkpoint
** 0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF 
**
** returns: refnum = number indicating the reference.
** \ingroup (CHKPT)
*/

int chkpt_rd_ref(void)
{
  int refnum;
  char *keyword;
  keyword = chkpt_build_keyword("Reference");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &refnum, sizeof(int) );

  free(keyword);
  return refnum;
}


/*!
** void chkpt_wt_ref(int)  
** Writes out the reference type from the flag in checkpoint
** 0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF 
**
** \param refnum = number indicating the reference.
** \ingroup (CHKPT)
*/

void chkpt_wt_ref(int refnum)
{
  char *keyword;
  keyword = chkpt_build_keyword("Reference");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &refnum, sizeof(int));

  free(keyword);
}
