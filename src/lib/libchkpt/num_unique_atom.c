/*!
  \file num_unique_atom.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_num_unique_atom()  
** Reads in the number of symmetry unique atoms.
**
** returns: nunique = number of symmetry unique atoms.
** \ingroup (CHKPT)
*/

int chkpt_rd_num_unique_atom(void)
{
  int nunique;

  psio_read_entry(PSIF_CHKPT, "::Num. unique atoms", (char *) &nunique, 
                  sizeof(int));
  return nunique;
}


/*!
** void chkpt_wt_num_unique_atom(int)  
** Writes out the number of symmetry unique atoms.
**
** \param nunique = number of symmetry unique atoms.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_num_unique_atom(int nunique)
{
  psio_write_entry(PSIF_CHKPT, "::Num. unique atoms", (char *) &nunique, 
                   sizeof(int));
}
