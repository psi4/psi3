/*!
  \file nsymhf.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nsymhf()  
** Reads in the total number of irreps in the point group 
** in which the molecule is being considered which 
** have non-zero number of basis functions.
**
** returns: int nirreps   total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
*/

int chkpt_rd_nsymhf(void)
{
  int nsymhf;

  psio_read_entry(PSIF_CHKPT, "::Num. HF irreps", (char *) &nsymhf, sizeof(int));
  return nsymhf;
}

/*!
** void chkpt_wt_nsymhf(int)  
** Writes out the total number of irreps in the point group 
** in which the molecule is being considered which 
** have non-zero number of basis functions.
**
** arguments: 
**  \param int nirreps   total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
*/

void chkpt_wt_nsymhf(int nsymhf)
{
  psio_write_entry(PSIF_CHKPT, "::Num. HF irreps", (char *) &nsymhf, sizeof(int));
}
