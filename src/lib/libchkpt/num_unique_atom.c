/*!
  \file num_unique_atom.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_num_unique_atom()  
** Reads in the number of symmetry unique atoms.
**
** returns: int nunique   number of symmetry unique atoms.
*/

int chkpt_rd_num_unique_atom(void)
{
  int nunique;

  psio_read_entry(PSIF_CHKPT, "::Num. unique atoms", (char *) &nunique, sizeof(int));
  return nunique;
}

/*!
** void chkpt_wt_num_unique_atom(int)  
** Writes out the number of symmetry unique atoms.
**
** arguments: 
**  \param int nunique   number of symmetry unique atoms.
**
** returns: none
*/

void chkpt_wt_num_unique_atom(int nunique)
{
  psio_write_entry(PSIF_CHKPT, "::Num. unique atoms", (char *) &nunique, sizeof(int));
}
