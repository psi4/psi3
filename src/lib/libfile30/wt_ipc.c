#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_ipc():	Writes a matrix of nshell*nirreps integers 
**			that contains SO transformation information.
**  FORMAT: ipc is stored in packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :                                 
**   | isc[3] | isc[2] | isc[1] | isc[0] |                 
** leftmost byte              rightmost byte
**			
**
**  arguments: int **ipc - integer matrix nshell by nirreps
**
**  returns: nothing	
*/


void file30_wt_ipc(int **ipc)
{
  int i;
  int *tmp;
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
  pack_4int(ipc, tmp, nshell, nirreps);
  
  wwritw(info30_.filenum, (char *) tmp, (int) nshell*npsym*sizeof(int),
	 ipc_ptr, &junk);

  free(tmp);
  return ;
}
