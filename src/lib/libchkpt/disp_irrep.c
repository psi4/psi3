/*!
  \file disp_irrep.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_disp_irrep()  
** Reads in the irrep of the current displaced geometry assuming
** Cotton ordering of irreps - to be used by input to determine
** docc and socc
**
** returns: disp_irrep = irrep of current displaced geometry
** \ingroup (CHKPT)
*/

int chkpt_rd_disp_irrep(void)
{
  int h;
  char *keyword;
  keyword = chkpt_build_keyword("Current Displacement Irrep");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &h, sizeof(int));

  free(keyword);
  return h;
}


/*!
** void chkpt_wt_disp_irrep(int)  
** Writes the irrep of the current displaced geometry assuming
** Cotton ordering of irreps - to be used by input to determine
** docc and socc
**
** \param disp_irrep = irrep of current displaced geometry
** \ingroup (CHKPT)
*/

void chkpt_wt_disp_irrep(int disp_irrep)
{
  char *keyword;
  keyword = chkpt_build_keyword("Current Displacement Irrep");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &disp_irrep, sizeof(int));

  free(keyword);
}
