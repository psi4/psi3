/*!
  \file rd_snuc.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_snuc(): Reads in array of the nuclei numbers shells belong to.
**
**  takes no arguments.
**
**  returns: int *snuc  an array of the nuclei numbers to which shells 
**                         belong to.
*/


int *file30_rd_snuc(void)
{
  int *snuc;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR snuc_ptr;

  nshell = file30_rd_nshell();
  snuc_ptr = (PSI_FPTR) (info30_.mpoint[7] - 1)*sizeof(int);

  snuc = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) snuc, (int) nshell*sizeof(int),
	 snuc_ptr, &junk);

  return snuc;
}
