/*!
  \file disp.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_disp():  Reads in the current geometry displacement number.
**
**   takes no arguments.
**
** not used by OPTKING; used by anybody else ???
**
**   returns: int disp   the current geometry displacement number
** \ingroup (CHKPT)
*/
int chkpt_rd_disp(void)
{
  int disp;
  char *keyword;
  keyword = chkpt_build_keyword("Current displacement");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &disp,
    sizeof(int));

  free(keyword);
  return disp;
}


/*!
** chkpt_wt_disp():  Writes out the current geometry displacement number.
**
**  arguments: 
**   \param int disp   the current geometry displacement number
**
** returns: none
** \ingroup (CHKPT)
*/
void chkpt_wt_disp(int disp)
{
  char *keyword;
  keyword = chkpt_build_keyword("Current displacement");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &disp,
    sizeof(int));

  free(keyword);
}
