/*!
  \file rd_sopi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_sopi():	
**
**  takes no arguments.
**
**  returns: int *sopi Read in an array nirreps long of the numbers of 
**			 symmetry orbitals in symm. blocks
*/


int *file30_rd_sopi(void)
{
  int *sopi;
  int nirreps;
  PSI_FPTR junk;
  PSI_FPTR sopi_ptr;

  nirreps = file30_rd_nirreps();
  sopi_ptr = (PSI_FPTR) (info30_.mpoint[36] - 1)*sizeof(int);

  sopi = init_int_array(nirreps);

  wreadw(info30_.filenum, (char *) sopi, (int) nirreps*sizeof(int),
	 sopi_ptr, &junk);

  return sopi;
}
