/*!
  \file num_unique_shell.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_num_unique_shell()  
** Reads in the number of symmetry unique shells. 
**
** returns: int nunique   number of symmetry unique shells.
*/


int chkpt_rd_num_unique_shell(void)
{
  int nunique;

  psio_read_entry(PSIF_CHKPT, "::Num. unique shells", (char *) &nunique, sizeof(int));
  return nunique;
}

/*!
** void chkpt_wt_num_unique_shell(int)  
** Writes out the number of symmetry unique shells. 
**
** arguments: 
**  \param int nunique   number of symmetry unique shells.
**
** returns: none
*/

void chkpt_wt_num_unique_shell(int nunique)
{
  psio_write_entry(PSIF_CHKPT, "::Num. unique shells", (char *) &nunique, sizeof(int));
}
