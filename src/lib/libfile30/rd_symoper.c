/*!
  \file rd_symoper.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_symoper():	Read in the mapping array between "canonical" ordering of symmetry
**                      operations in the point group and the one defined in symmetry.h
**
**  takes no arguments.
**
**  returns: int *symoper Array nirrep long
*/


int *file30_rd_symoper(void)
{
  int *symoper;
  int nirreps;
  PSI_FPTR symoper_ptr;

  nirreps = file30_rd_nirreps();
  symoper_ptr = (PSI_FPTR) (info30_.mpoint[44] - 1)*sizeof(int);

  symoper = init_int_array(nirreps);

  wreadw(info30_.filenum, (char *) symoper, (int) nirreps*sizeof(int),
	 symoper_ptr, &symoper_ptr);

  return symoper;
}
