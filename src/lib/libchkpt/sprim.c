/*!
  \file sprim.c
*/

#include <stdio.h>
#include <libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sprim():	Reads in array of the numbers of first primitives 
**			from the shells.
**
**  takes no arguments.
**
**  returns: int *sprim an array of the numbers of first primitives
**			from the shells.
*/


int *chkpt_rd_sprim(void)
{
  int *sprim;
  int nshell;

  nshell = chkpt_rd_nshell();

  sprim = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::First primitive per shell", (char *) sprim, 
		  nshell*sizeof(int));

  return sprim;
}

/*!
** chkpt_wt_sprim():	Writes out an array of the numbers of first primitives 
**			from the shells.
**
**  arguments:
**  \param int *: sprim: an array of the numbers of first primitives
**			from the shells.
**
**  returns: none
*/


void chkpt_wt_sprim(int *sprim)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::First primitive per shell", (char *) sprim, 
		   nshell*sizeof(int));
}
