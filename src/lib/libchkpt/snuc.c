/*!
  \file snuc.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_snuc(): Reads in array of the nuclei numbers shells belong to.
**
**  takes no arguments.
**
**  returns: snuc = an array of the nuclei numbers to which shells 
**                  belong to.
**
** \ingroup (CHKPT)
*/

int *chkpt_rd_snuc(void)
{
  int *snuc;
  int nshell;

  nshell = chkpt_rd_nshell();

  snuc = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::Shell nucleus", (char *) snuc, 
                  nshell*sizeof(int));

  return snuc;
}


/*!
** chkpt_wt_snuc(): Writes out array of the nuclei numbers shells belong to.
**
**  \param snuc = an array of the nuclei numbers to which shells belong to
**
**  returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_snuc(int *snuc)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::Shell nucleus", (char *) snuc, 
                   nshell*sizeof(int));
}
