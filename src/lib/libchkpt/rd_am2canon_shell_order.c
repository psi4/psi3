/*!
  \file rd_am2canon_shell_order.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_am2canon_shell_order() 
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** returns: int *am2can_shell_order
*/


int *chkpt_rd_am2canon_shell_order(void)
{
  int *am2can_sh_ord, nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  am2can_sh_ord = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Am2can sh ord", (char *) am2can_sh_ord, nshell*sizeof(int),
            next, &next); 

  return am2can_sh_ord;
}
