/*!
  \file nshell.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nshell() 
** Reads in the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells. 
** Shells of all atoms are counted (compare nprim).
**
** returns: nshell = total number of shells.
** \ingroup (CHKPT)
*/

int chkpt_rd_nshell(void)
{
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("Num. shells");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nshell, sizeof(int));

  free(keyword);
  return nshell;
}


/*!
** void chkpt_wt_nshell(int) 
** Writes out the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells. 
** Shells of all atoms are counted (compare nprim).
**
** \param nshell = total number of shells.
**
** returns:none 
**
** \ingroup (CHKPT)
*/

void chkpt_wt_nshell(int nshell)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. shells");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nshell, sizeof(int));

  free(keyword);
}
