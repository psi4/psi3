/*!
  \file rd_sprim.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_sprim()
** Reads in array of the numbers of first primitives 
** from the shells.
**
** returns: int *sprim  an array of the numbers of first primitives
**			from the shells.
*/


int *chkpt_rd_sprim(void)
{
  int *sprim;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  sprim = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Sprim", (char *) sprim, nshell*sizeof(int),
            next, &next);

  return sprim;
}
