/*!
  \file rd_nprim.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_nprim()  
** Reads in the total number of primitive Gaussian functions 
** (only primitives of symmetry independent atoms are taken into account!).
**
**   returns: int nprim   total number of primitive Gaussian functions.
*/


int chkpt_rd_nprim(void)
{
  int nprim;

  psio_read_entry(PSIF_CHKPT, "::Num prim", (char *) &nprim, 
                  sizeof(int) );
  return nprim;
}
