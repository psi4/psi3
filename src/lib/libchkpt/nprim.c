/*!
  \file nprim.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
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

  psio_read_entry(PSIF_CHKPT, "::Num. primitives", (char *) &nprim, 
                  sizeof(int));
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
  psio_write_entry(PSIF_CHKPT, "::Num. primitives", (char *) &nprim, 
                   sizeof(int));
}
