/*!
  \file rd_symoper.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_symoper()
** Read in the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
**  returns: int *symoper    Array nirrep long
*/


int *chkpt_rd_symoper(void)
{
  int *symoper;
  int nirreps;
  psio_address next;
  next = PSIO_ZERO;

  nirreps = chkpt_rd_nirreps();
  symoper = init_int_array(nirreps);

  psio_read(PSIF_CHKPT, "::Symoper", (char *) symoper, nirreps*sizeof(int),
            next, &next);

  return symoper;
}
