/*!
  \file nallatom.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_nallatom()  
**
** Reads number of all atoms (including dummy atoms)
**
** takes no arguments.
**
** returns: nallatom = number of all atoms (including dummies).
** \ingroup (CHKPT)
*/

int chkpt_rd_nallatom(void)
{
  int num_allatoms;
  char *keyword;
  keyword = chkpt_build_keyword("Num. all atoms");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms, sizeof(int));

  free(keyword);
  return num_allatoms;
}


/*!
** chkpt_wt_nallatom()
**
** Writes the number of all atoms (including dummy atoms)
**
** \param nallatom  = number of all atoms (including dummies).
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_nallatom(int num_allatoms)
{
  char *keyword;
  keyword = chkpt_build_keyword("Num. all atoms");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms, sizeof(int));

  free(keyword);
}
