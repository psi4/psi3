/*!
  \file snumg.c
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_snumg():	Reads in array of the numbers of the primitive Gaussians
**			in shells.
**
**  takes no arguments.
**
**  returns: int *snumg	Reads in array of the numbers of the primitive Gaussians
**			in shells
*/


int *chkpt_rd_snumg(void)
{
  int *snumg;
  int nshell;

  nshell = chkpt_rd_nshell();
  snumg = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::Primitives per shell", (char *) snumg, nshell*sizeof(int));

  return snumg;
}

/*!
** chkpt_wt_snumg():	Writes out array of the numbers of the primitive Gaussians
**			in shells.
**
**  arguments:
**  \param int *snumg:	 array of the numbers of the primitive Gaussians
**			in shells
*/


void chkpt_wt_snumg(int *snumg)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::Primitives per shell", (char *) snumg, nshell*sizeof(int));
}
