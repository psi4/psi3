/*!
  \file rd_snumg.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_snumg():	Reads in array of the numbers of the primitive Gaussians
**			in shells.
**
**  takes no arguments.
**
**  returns: int *snumg	Reads in array of the numbers of the primitive Gaussians
**			in shells
*/


int *file30_rd_snumg(void)
{
  int *snumg;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR snumg_ptr;

  nshell = file30_rd_nshell();
  snumg_ptr = (PSI_FPTR) (info30_.mpoint[9] - 1)*sizeof(int);

  snumg = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) snumg, (int) nshell*sizeof(int),
	 snumg_ptr, &junk);

  return snumg;
}
