/*!
  \file rd_us2c.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_us2s()
** Read in a mapping array betwen unique shell and 
** full shell lists
**
**  returns: int *us2s   Read in an array num_unique_shell
*/


int *chkpt_rd_us2s(void)
{
  int *us2s;
  int num_unique_shells;
  psio_address next;
  next = PSIO_ZERO;

  num_unique_shells = chkpt_rd_num_unique_shell();
  us2s = init_int_array(num_unique_shells);

  psio_read(PSIF_CHKPT, "::Us2s", (char *) us2s, num_unique_shells*sizeof(int),
            next, &next);

  return us2s;
}
