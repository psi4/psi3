/*!
  \file sprim.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sprim(): Reads in array of the numbers of first primitives 
**                   from the shells.
**
**  takes no arguments.
**
**  returns: sprim = an array of the numbers of first primitives
**           from the shells.
**
** \ingroup (CHKPT)
*/

int *chkpt_rd_sprim(void)
{
  int *sprim;
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("First primitive per shell");

  nshell = chkpt_rd_nshell();

  sprim = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

  free(keyword);
  return sprim;
}


/*!
** chkpt_wt_sprim():	Writes out an array of the numbers of first primitives 
**			from the shells.
**
**  \param sprim = an array of the numbers of first primitives
**                 from the shells.
**
**  returns: none
** 
** \ingroup (CHKPT)
*/

void chkpt_wt_sprim(int *sprim)
{
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("First primitive per shell");

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

  free(keyword);
}
