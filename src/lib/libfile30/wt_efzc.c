#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_efzc():  Writes the frozen-core energy to file30.
**
**    arguments: double efzc  the frozen-core energy.
**
*/

void file30_wt_efzc(double efzc)
{
  int natom;
  PSI_FPTR junk;
  PSI_FPTR ptr;

  natom = file30_rd_natom();
  ptr = (PSI_FPTR) 
        ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+4*sizeof(double));

  wwritw(info30_.filenum, (char *) &efzc, (int) sizeof(double),
	 ptr, &junk);

}
