/*!
  \file nao.c
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

  psio_read_entry(PSIF_CHKPT, "::Num. AO", (char *) &nao, sizeof(int));
  return nao;
}

/*!
**  void chkpt_wt_nao(int)
**  Writes out the total number of atomic orbitals.
**
**  arguments: 
**   \param int nao   total number of atomic orbitals.
**
**  returns: none
*/

void chkpt_wt_nao(int nao)
{
  psio_write_entry(PSIF_CHKPT, "::Num. AO", (char *) &nao, sizeof(int));
}
