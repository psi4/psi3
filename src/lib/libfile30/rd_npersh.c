#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_npersh(): 
**
**  takes no arguments.
**
**  returns: int npersh[MAXANGMOM]
*/


int *file30_rd_npersh(void)
{
  int *npersh;
  PSI_FPTR junk;
  PSI_FPTR npersh_ptr;

  npersh_ptr = (PSI_FPTR) (info30_.mpoint[35] - 1)*sizeof(int);

  npersh = init_int_array(MAXANGMOM);

  wreadw(info30_.filenum, (char *) npersh, (int) MAXANGMOM*sizeof(int),
	 npersh_ptr, &junk);

  return npersh;
}
