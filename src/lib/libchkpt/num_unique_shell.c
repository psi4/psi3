/*!
  \file num_unique_shell.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. unique shells");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

  free(keyword);
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. unique shells");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

  free(keyword);
}
