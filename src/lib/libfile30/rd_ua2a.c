/*!
  \file rd_ua2a.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_ua2a(): Read in a mapping array from the symmetry-unique atom
**                   list to the full atom list
**
**  takes no arguments.
**
**  returns: int *ua2a Read in an array num_unique_atom long
*/


int *file30_rd_ua2a(void)
{
  int *ua2a;
  int num_unique_atoms;
  PSI_FPTR ua2a_ptr;

  num_unique_atoms = file30_rd_num_unique_atom();
  ua2a_ptr = (PSI_FPTR) (info30_.mpoint[43] - 1)*sizeof(int);

  ua2a = init_int_array(num_unique_atoms);

  wreadw(info30_.filenum, (char *) ua2a, (int) num_unique_atoms*sizeof(int),
	 ua2a_ptr, &ua2a_ptr);

  return ua2a;
}
