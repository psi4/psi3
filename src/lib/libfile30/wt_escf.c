#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_escf(): Writes the scf energy to file30.
**
**   arguments: double escf  the scf energy.
**
**   returns nothing.
*/

void file30_wt_escf(double escf)
{
  int  natom;
  PSI_FPTR junk;
  PSI_FPTR escf_ptr;

  natom = file30_rd_natom();
  escf_ptr = (PSI_FPTR) 
           ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+sizeof(double));

  wwritw(info30_.filenum, (char *) &escf, (int) sizeof(double),
	 escf_ptr, &junk);

}
