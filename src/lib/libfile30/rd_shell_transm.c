#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_shell_transm():	Read in a matrix of nshell*nirreps integers 
**			        that contains symmetry information.
**
**  takes no arguments.
**
**  returns: int **shell_transm	
*/


int **file30_rd_shell_transm(void)
{
  int i, shell, irrep;
  int **shell_transm;
  int nshell;
  int nirreps;
  PSI_FPTR shell_transm_ptr;

  nshell = file30_rd_nshell();
  nirreps = file30_rd_nirreps();
  shell_transm_ptr = (PSI_FPTR) (info30_.mpoint[26] - 1)*sizeof(int);

  shell_transm = init_int_matrix(nshell,nirreps);
  for(shell=0;shell<nshell;shell++)
    wreadw(info30_.filenum, (char *) shell_transm[shell], nirreps*sizeof(int),
	   shell_transm_ptr, &shell_transm_ptr);

  return shell_transm;
}
