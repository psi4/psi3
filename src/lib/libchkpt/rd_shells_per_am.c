/*!
  \file rd_shell_per_am.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_shells_per_am() 
** Reads in the numbers of shells of each angular momentum.
**
**  returns: int *shells_per_am
*/


int *chkpt_rd_shells_per_am(void)
{
  int *shells_per_am;
  int max_am;
  psio_address next;
  next = PSIO_ZERO;

  max_am = chkpt_rd_max_am();
  shells_per_am = init_int_array(max_am+1);

  psio_read(PSIF_CHKPT, "::Shells per am", (char *) shells_per_am,
            (max_am+1)*sizeof(int), next, &next);

  return shells_per_am; 
}
