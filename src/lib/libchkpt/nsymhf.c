/*!
  \file nsymhf.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nsymhf()  
** Reads in the total number of irreps in the point group 
** in which the molecule is being considered which 
** have non-zero number of basis functions.
**
** returns: nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup (CHKPT)
*/

int chkpt_rd_nsymhf(void)
{
  int nsymhf;
  char *keyword;
  keyword = chkpt_build_keyword("Num. HF irreps");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

  free(keyword);
  return nsymhf;
}


/*!
** void chkpt_wt_nsymhf(int)  
** Writes out the total number of irreps in the point group 
** in which the molecule is being considered which 
** have non-zero number of basis functions.
**
** \param nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup (CHKPT)
*/

void chkpt_wt_nsymhf(int nsymhf)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. HF irreps");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

  free(keyword);
}
