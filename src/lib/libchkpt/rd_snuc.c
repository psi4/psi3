/*!
  \file rd_snuc.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_snuc() 
** Reads in array of the nuclei numbers shells belong to.
**
**  returns: int *snuc  an array of the nuclei numbers to which shells 
**                      belong to.
*/


int *chkpt_rd_snuc(void)
{
  int *snuc;
  int nshell;
  psio_address next;
  next = PSIO_ZERO;

  nshell = chkpt_rd_nshell();
  snuc = init_int_array(nshell);

  psio_read(PSIF_CHKPT, "::Snuc", (char *) snuc, nshell*sizeof(int),
            next, &next);

  return snuc;
}
