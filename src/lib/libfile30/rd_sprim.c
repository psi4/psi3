#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_sprim():	Reads in array of the numbers of first primitives 
**			from the shells.
**
**  takes no arguments.
**
**  returns: int *sprim an array of the numbers of first primitives
**			from the shells.
*/


int *file30_rd_sprim(void)
{
  int *sprim;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR sprim_ptr;

  nshell = file30_rd_nshell();
  sprim_ptr = (PSI_FPTR) (info30_.mpoint[6] - 1)*sizeof(int);

  sprim = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) sprim, (int) nshell*sizeof(int),
	 sprim_ptr, &junk);

  return sprim;
}
