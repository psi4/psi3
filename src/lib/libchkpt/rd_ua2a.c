/*!
  \file rd_ua2a.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_ua2a()
** Read in a mapping array from the symmetry-unique atom
** list to the full atom list
**
**  returns: int *ua2a    Read in an array num_unique_atom long
*/


int *chkpt_rd_ua2a(void)
{
  int *ua2a;
  int num_unique_atoms;
  psio_address next;
  next = PSIO_ZERO;

  num_unique_atoms = chkpt_rd_num_unique_atom();
  ua2a = init_int_array(num_unique_atoms);

  psio_read(PSIF_CHKPT, "::Ua2a", (char *) ua2a, num_unique_atoms*sizeof(int),
            next, &next);

  return ua2a;
}
