#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_sloc():	Read in an array of the numbers of the first AO 
**			from the shells.
**
**  takes no arguments.
**
**  returns: int *sloc	Read in an array nshell long of the numbers of 
**			the first AOs from the shells.
*/


int *file30_rd_sloc(void)
{
  int *sloc;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR sloc_ptr;

  nshell = file30_rd_nshell();
  sloc_ptr = (PSI_FPTR) (info30_.mpoint[10] - 1)*sizeof(int);

  sloc = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) sloc, (int) nshell*sizeof(int),
	 sloc_ptr, &junk);

  return sloc;
}
