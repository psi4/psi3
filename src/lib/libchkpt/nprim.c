/*!
  \file nprim.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nprim()  
** Reads in the total number of primitive Gaussian functions 
** (only primitives of symmetry independent atoms are taken into account!).
**
** returns: nprim = total number of primitive Gaussian functions.
** \ingroup (CHKPT)
*/

int chkpt_rd_nprim(void)
{
  int nprim;
  char *keyword;
  keyword = chkpt_build_keyword("Num. primitives");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

  free(keyword);
  return nprim;
}


/*!
** void chkpt_wt_nprim(int)  
** Writes out the total number of primitive Gaussian functions 
** (only primitives of symmetry independent atoms are taken into account!).
**
** \param nprim = total number of primitive Gaussian functions.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_nprim(int nprim)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. primitives");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

  free(keyword);
}
