/*!
  \file rd_num_unique_atom.c
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

  psio_read_entry(PSIF_CHKPT, "::Num unique atom", (char *) &nunique, 
                  sizeof(int) );
  return nunique;
}
