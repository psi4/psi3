#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_irrstr():	Read in an array of irreps' characters packed in ints
**
**  takes no arguments.
**
**  returns: int *irrstr Array nirrep long
*/


int *file30_rd_irrstr(void)
{
  int *irrstr;
  int nirreps;
  PSI_FPTR irrstr_ptr;

  nirreps = file30_rd_nirreps();
  irrstr_ptr = (PSI_FPTR) (info30_.mpoint[43] - 1)*sizeof(int);

  irrstr = init_int_array(nirreps);

  wreadw(info30_.filenum, (char *) irrstr, (int) nirreps*sizeof(int),
	 irrstr_ptr, &irrstr_ptr);

  return irrstr;
}
