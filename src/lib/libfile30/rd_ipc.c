#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_ipc():	Reads in a matrix of SO transformation information
**			( see libfile30.h for details)
**
**  FORMAT: ipc is stored in packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :
**   | ipc[3] | ipc[2] | ipc[1] | ipc[0] |
** leftmost byte              rightmost byte
**
**  takes no arguments.
**
**  returns: int **ipc	Returns a matrix of nshell*nirreps integers
**                      ( see libfile30.h for details).
*/


int **file30_rd_ipc(void)
{
  int i,j,count;
  int **ipc, **packed_ipc, *tmp;
  int nshell;
  int nirreps;
  int npsym;
  PSI_FPTR junk;
  PSI_FPTR ipc_ptr;

  nshell = file30_rd_nshell();
  nirreps = file30_rd_nirreps();
  ipc_ptr = (PSI_FPTR) (info30_.mpoint[14] - 1)*sizeof(int);

  npsym = (nirreps == 8) ? 2 : 1;
  tmp = init_int_array(nshell*npsym);

  wreadw(info30_.filenum, (char *) tmp, (int) nshell*npsym*sizeof(int),
	 ipc_ptr, &junk);

  ipc = init_int_matrix(nshell,nirreps);
  unpack_4int(tmp, ipc, nshell, nirreps);

  free(tmp);
  return ipc;
}
