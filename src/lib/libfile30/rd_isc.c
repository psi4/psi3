#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_isc():	Read in a matrix of nshell*nirreps integers 
**			that contains symmetry information.
**  FORMAT: isc is stored in packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :                                 
**   | isc[3] | isc[2] | isc[1] | isc[0] |                 
** leftmost byte              rightmost byte
**			
**
**  takes no arguments.
**
**  returns: int **isc	
*/


int **file30_rd_isc(void)
{
  int i;
  int **isc, *tmp;
  int nshell;
  int nirreps;
  int npsym;
  PSI_FPTR junk;
  PSI_FPTR isc_ptr;

  nshell = file30_rd_nshell();
  nirreps = file30_rd_nirreps();
  isc_ptr = (PSI_FPTR) (info30_.mpoint[13] - 1)*sizeof(int);

  npsym = (nirreps == 8) ? 2 : 1;
  tmp = init_int_array(nshell*npsym);

  wreadw(info30_.filenum, (char *) tmp, (int) nshell*npsym*sizeof(int),
	 isc_ptr, &junk);

  isc = init_int_matrix(nshell,nirreps);
  unpack_4int(tmp, isc, nshell, nirreps);
  
  free(tmp);
  return isc;
}
