/*!
  \file rd_num_unique_shell.c
*/

#include "chkpt.h"
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

  psio_read_entry(PSIF_CHKPT, "::Num unique shell", (char *) &nunique, 
                  sizeof(int) );
  return nunique;
}
