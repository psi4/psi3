#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_smax():	Reads in the last index in the offset array for 
** 	that type of the shell. That array used for manipulation with 
**	indices in integral calculations.
**	For s-type shells this index is 1(there exists only one possible
**	s-type cartesian primitive Gaussian), for p-type - 4 (there are 3 
**	possible p-type primitives), for d-type - 10, for f-type - 20 ... 
**	This index exists thanks to the limit on the possible angular momentum
**	of the basis functions. Hopefully, it will be changed soon.
**	Note : smax(l) = smax(l-1) + (l+2)*(l+1)/2
**
**  takes no arguments.
**
**  returns: int *smax   The indices are returned as an array of integers.
*/


int *file30_rd_smax(void)
{
  int *smax;
  int nshell = 0;
  PSI_FPTR junk;
  PSI_FPTR smax_ptr;

  nshell = file30_rd_nshell();
  smax_ptr = (PSI_FPTR) (info30_.mpoint[12] - 1)*sizeof(int);

  smax = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) smax, (int) nshell*sizeof(int),
	 smax_ptr, &junk);

  return smax;
}
