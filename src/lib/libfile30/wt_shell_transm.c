#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_shell_transm():	Write out a matrix of nshell*nirreps integers 
**			        that contains shell orbits.
**  FORMAT: shell_transm is stored in packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :                                 
**   | shell_transm[3] | shell_transm[2] | shell_transm[1] | shell_transm[0] |
**     leftmost byte                                         rightmost byte
**			
**
**  arguments: int **shell_transm
**
**  returns nothing
*/


void file30_wt_shell_transm(int **shell_transm)
{
  int i;
  int *tmp;
  int nshell;
  int nirreps;
  int npsym;
  PSI_FPTR junk;
  PSI_FPTR shell_transm_ptr;

  nshell = file30_rd_nshell();
  nirreps = file30_rd_nirreps();
  shell_transm_ptr = (PSI_FPTR) (info30_.mpoint[14] - 1)*sizeof(int);

  npsym = (nirreps == 8) ? 2 : 1;
  tmp = init_int_array(nshell*npsym);
  pack_4int(shell_transm, tmp, nshell, nirreps);
  
  wwritw(info30_.filenum, (char *) tmp, (int) nshell*npsym*sizeof(int),
	 shell_transm_ptr, &junk);

  free(tmp);
  return ;
}
