#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_smin():	Reads in the first index in the offset array for 
** 	that type of the shell. That array used for manipulation with 
**	indices in integral calculations.
**	For s-type shells this index is 1, for p-type - 2 (there exists 
**	only one possible s-type cartesian primitive Gaussian), for d-type - 5
**	(there are 3 possible p-type primitives), for f-type - 11 ... This
**	index exists thanks to the limit on the possible angular momentum
**	of the basis functions. Hopefully, it will be changed soon.
**	Note : smin(l) = smin(l-1) + (l+1)*l/2
**
**  takes no arguments.
**
**  returns: int *smin   The indices are returned as an array of integers.
*/


int *file30_rd_smin(void)
{
  int *smin;
  int nshell = 0;
  PSI_FPTR junk;
  PSI_FPTR smin_ptr;

  nshell = file30_rd_nshell();
  smin_ptr = (PSI_FPTR) (info30_.mpoint[11] - 1)*sizeof(int);

  smin = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) smin, (int) nshell*sizeof(int),
	 smin_ptr, &junk);

  return smin;
}
