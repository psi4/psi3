/*!
  \file rd_us2c.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_us2s(): Read in a mapping array betwen unique shell and 
**		     full shell lists
**
**  takes no arguments.
**
**  returns: int *us2s Read in an array num_unique_shell
*/


int *file30_rd_us2s(void)
{
  int *us2s;
  int num_unique_shells;
  PSI_FPTR junk;
  PSI_FPTR us2s_ptr;

  num_unique_shells = file30_rd_num_unique_shell();
  us2s_ptr = (PSI_FPTR) (info30_.mpoint[39] - 1)*sizeof(int);

  us2s = init_int_array(num_unique_shells);

  wreadw(info30_.filenum, (char *) us2s, (int) num_unique_shells*sizeof(int),
	 us2s_ptr, &junk);

  return us2s;
}
