/*!
  \file rd_nao.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
**  int chkpt_rd_nao()
**  Reads in the total number of atomic orbitals.
**
**  returns: int nao   total number of atomic orbitals.
*/


int chkpt_rd_nao(void)
{
  int nao;

  psio_read_entry(PSIF_CHKPT, "::Num ao", (char *) &nao, 
                  sizeof(int) );
  return nao;
}
