/*!
  \file rd_snumg.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_snumg()	
** Reads in array of the numbers of the primitive Gaussians
** in shells.
**
** returns: int *snumg	  Reads in array of the numbers of the primitive Gaussians
*                         in shells
*/


int *chkpt_rd_snumg(void)
{
  int *snumg;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  snumg = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Snumg", (char *) snumg, nshell*sizeof(int),
            next, &next);

  return snumg;
}
