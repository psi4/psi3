/*!
  \file rd_nshell.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nshell() 
** Reads in the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells. 
** Shells of all atoms are counted (compare nprim).
**
**   returns: int nshell   total number of shells.
*/


int chkpt_rd_nshell(void)
{
  int nshell;

  psio_read_entry(PSIF_CHKPT, "::Num shells", (char *) &nshell, 
                  sizeof(int) );
  return nshell;
}
