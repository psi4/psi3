/*!
  \file sloc.c
*/

#include <stdio.h>
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
**  returns: int *sloc	An array nshell long of the numbers of 
**			the first AOs from the shells.
*/


int *chkpt_rd_sloc(void)
{
  int *sloc;
  int nshell;

  nshell = chkpt_rd_nshell();

  sloc = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::First AO per shell", (char *) sloc, nshell*sizeof(int));

  return sloc;
}

/*!
** chkpt_wt_sloc():	Writes out an array of the numbers of the first AO 
**			from the shells.
**
**  arguments: 
**  \param int *sloc:	An array nshell long of the numbers of 
**			the first AOs from the shells.
**  returns: none
*/


void chkpt_wt_sloc(int *sloc)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::First AO per shell", (char *) sloc, nshell*sizeof(int));
}
