/*!
  \file num_unique_shell.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_num_unique_shell()  
** Reads in the number of symmetry unique shells. 
**
** returns: nunique = number of symmetry unique shells.
** \ingroup (CHKPT)
*/

int chkpt_rd_num_unique_shell(void)
{
  int nunique;

  psio_read_entry(PSIF_CHKPT, "::Num. unique shells", (char *) &nunique, 
                  sizeof(int));
  return nunique;
}


/*!
** void chkpt_wt_num_unique_shell(int)  
** Writes out the number of symmetry unique shells. 
**
** \param nunique = number of symmetry unique shells.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_num_unique_shell(int nunique)
{
  psio_write_entry(PSIF_CHKPT, "::Num. unique shells", (char *) &nunique, 
                   sizeof(int));
}
