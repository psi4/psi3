/*!
  \file nirreps.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nirreps()  
** Reads in the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** returns: nirreps = total number of irreducible representations.
** \ingroup (CHKPT)
*/

int chkpt_rd_nirreps(void)
{
  int nirreps;

  psio_read_entry(PSIF_CHKPT, "::Num. irreps", (char *) &nirreps, sizeof(int));
  return nirreps;
}


/*!
** void chkpt_wt_nirreps(int)  
** Writes out the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** \param nirreps = total number of irreducible representations.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_nirreps(int nirreps)
{
  psio_write_entry(PSIF_CHKPT, "::Num. irreps", (char *) &nirreps, sizeof(int));
}
