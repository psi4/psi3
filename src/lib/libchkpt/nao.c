/*!
  \file nao.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
**  int chkpt_rd_nao()
**  Reads in the total number of atomic orbitals.
**
**  returns: nao = total number of atomic orbitals.
**  \ingroup (CHKPT)
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
**   \param nao = total number of atomic orbitals.
**
**  returns: none
**  \ingroup (CHKPT)
*/

void chkpt_wt_nao(int nao)
{
  psio_write_entry(PSIF_CHKPT, "::Num. AO", (char *) &nao, sizeof(int));
}
