#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_am2canon_shell_order(): Reads in the mapping array from the am-ordered
**                                   to the canonical (in the order of appearance)
**                                   list of shells.
**
**  takes no arguments.
**
**  returns: int *am2canon_shell_order
*/


int *file30_rd_am2canon_shell_order(void)
{
  int *am2can_sh_ord;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR ptr;

  nshell = file30_rd_nshell();
  ptr = (PSI_FPTR) (info30_.mpoint[48] - 1)*sizeof(int);

  am2can_sh_ord = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) am2can_sh_ord, (int) nshell*sizeof(int),
	 ptr, &junk);

  return am2can_sh_ord;
}
