/*!
  \file sloc.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sloc():	Read in an array of the numbers of the first AO 
**			from the shells.
**
**  takes no arguments.
**
**  returns: sloc = An array nshell long of the numbers of 
**                  the first AOs from the shells.
**
** \ingroup (CHKPT)
*/

int *chkpt_rd_sloc(void)
{
  int *sloc;
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("First AO per shell");

  nshell = chkpt_rd_nshell();

  sloc = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

  free(keyword);
  return sloc;
}


/*!
** chkpt_wt_sloc():	
**  Writes out an array of the numbers of the first AO from the shells.
**
**  \param sloc = An array nshell long of the numbers of the first AOs 
**                from the shells.
**  returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_sloc(int *sloc)
{
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("First AO per shell");

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

  free(keyword);
}
