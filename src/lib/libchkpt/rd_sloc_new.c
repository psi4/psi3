/*!
  \file rd_sloc_new.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_sloc_new()	
** Read in an array of the numbers of the first basis 
** functions (not AOs as rd_sloc does)  from the shells.
**
** returns: int *sloc	Read in an array nshell long of the numbers of 
**			the first basis functions from the shells.
*/


int *chkpt_rd_sloc_new(void)
{
  int *sloc_new;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  sloc_new = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Sloc new", (char *) sloc_new, nshell*sizeof(int),
            next, &next);

  return sloc_new;
}
