/*!
  \file rd_nirreps.c
*/

#include "chkpt.h"
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nirreps()  
** Reads in the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** returns: int nirreps   total number of irreducible representations.
*/

int chkpt_rd_nirreps(void)
{
  int nirreps;

  psio_read_entry(PSIF_CHKPT, "::Num irreps", (char *) &nirreps, 
                  sizeof(int) );
  return nirreps;
}
