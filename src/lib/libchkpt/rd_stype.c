/*!
  \file rd_stype.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_stype()
** Reads in an array of the angular momentum numbers of 
** the shells.
**
** returns: int *stype	an array of the angular momentum numbers of the
**			shells.
*/


int *chkpt_rd_stype(void)
{
  int *stype;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  stype = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Stype", (char *) stype, nshell*sizeof(int),
            next, &next);

  return stype;
}
