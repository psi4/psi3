/*!
  \file nao.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. AO");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &nao, sizeof(int));

  free(keyword);
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
  char *keyword;
  keyword = chkpt_build_keyword("Num. AO");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &nao, sizeof(int));

  free(keyword);
}
