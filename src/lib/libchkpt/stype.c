/*!
  \file stype.c
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_stype(): 	Reads in an array of the angular momentum numbers of 
**			the shells.
**
**  takes no arguments.
**
**  returns: int *stype	an array of the angular momentum numbers of the
**			shells.
*/


int *chkpt_rd_stype(void)
{
  int *stype;
  int nshell;

  nshell = chkpt_rd_nshell();

  stype = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::Shell ang. mom.", (char *) stype, nshell*sizeof(int));

  return stype;
}

/*!
** chkpt_wt_stype(): 	Writes out an array of the angular momentum numbers of 
**			the shells.
**
**  arguments: 
**  \param int *stype:	an array of the angular momentum numbers of the
**			shells.
**  returns: none
*/


void chkpt_wt_stype(int *stype)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::Shell ang. mom.", (char *) stype, nshell*sizeof(int));
}
