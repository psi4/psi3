/*!
  \file rd_sloc.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_sloc()	
** Read in an array of the numbers of the first AO 
** from the shells.
**
** returns: int *sloc	Read in an array nshell long of the numbers of 
**			the first AOs from the shells.
*/


int *chkpt_rd_sloc(void)
{
  int *sloc;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  sloc = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Sloc", (char *) sloc, nshell*sizeof(int),
            next, &next);

  return sloc;
}
