/*!
  \file rd_stype.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_stype(): 	Reads in an array of the angular momentum numbers of 
**			the shells.
**
**  takes no arguments.
**
**  returns: int *stype	an array of the angular momentum numbers of the
**			shells.
*/


int *file30_rd_stype(void)
{
  int *stype;
  int nshell;
  PSI_FPTR junk;
  PSI_FPTR stype_ptr;

  nshell = file30_rd_nshell();
  stype_ptr = (PSI_FPTR) (info30_.mpoint[8] - 1)*sizeof(int);

  stype = init_int_array(nshell);

  wreadw(info30_.filenum, (char *) stype, (int) nshell*sizeof(int),
	 stype_ptr, &junk);

  return stype;
}
