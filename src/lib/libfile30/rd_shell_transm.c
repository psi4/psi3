#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_shell_transm():	Read in a matrix of nshell*nirreps integers 
**			        that contains symmetry information.
**  FORMAT: shell_transm is stored in packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :                                 
**   | shell_transm[3] | shell_transm[2] | shell_transm[1] | shell_transm[0] |
**     leftmost byte                                         rightmost byte
**			
**
**  takes no arguments.
**
**  returns: int **shell_transm	
*/


int **file30_rd_shell_transm(void)
{
  int i;
  int **shell_transm, *tmp;
  int nshell;
  int nirreps;
  int npsym;
  PSI_FPTR junk;
  PSI_FPTR shell_transm_ptr;

  nshell = file30_rd_nshell();
  nirreps = file30_rd_nirreps();
  shell_transm_ptr = (PSI_FPTR) (info30_.mpoint[26] - 1)*sizeof(int);

  npsym = (nirreps == 8) ? 2 : 1;
  tmp = init_int_array(nshell*npsym);

  wreadw(info30_.filenum, (char *) tmp, (int) nshell*npsym*sizeof(int),
	 shell_transm_ptr, &junk);

  shell_transm = init_int_matrix(nshell,nirreps);
  unpack_4int(tmp, shell_transm, nshell, nirreps);

  free(tmp);
  return shell_transm;
}
