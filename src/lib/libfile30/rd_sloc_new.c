/*!
  \file rd_sloc_new.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_sloc_new():	Read in an array of the numbers of the first basis 
**			functions (not AOs as rd_sloc does)  from the shells.
**
**  takes no arguments.
**
**  returns: int *sloc	Read in an array nshell long of the numbers of 
**			the first basis functions from the shells.
*/


int *file30_rd_sloc_new(void)
{
  int *sloc_new;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR sloc_new_ptr;

  nshell = file30_rd_nshell();
  sloc_new_ptr = (PSI_FPTR) (info30_.mpoint[42] - 1)*sizeof(int);

  sloc_new = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) sloc_new, (int) nshell*sizeof(int),
	 sloc_new_ptr, &junk);

  return sloc_new;
}
