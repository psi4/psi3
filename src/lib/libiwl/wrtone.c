#include <stdio.h>
#include <psio.h>
#include "iwl.h"

/*
** IWL_WRTONE()
**
** This function writes the one-electron integrals.
**
** Arguments:
**   itap       = tape to read ints from
**   ntri       = the size of the array (lower triangle)
**   onel_ints  = array to hold the one-electron integrals.
**   e_fzc      = frozen core energy
**
** David Sherrill, March 1995
*/
void iwl_wrtone(int itap, int ntri, double *onel_ints, double e_fzc)
{
  psio_open(itap,PSIO_OPEN_NEW);  /* We assume that we can overwrite */
  psio_write_entry(itap, IWL_KEY_EFZC, (char *) &e_fzc, sizeof(double));
  psio_write_entry(itap, IWL_KEY_ONEL, (char *) onel_ints, ntri*sizeof(double));
  psio_close(itap,1);
}


