/*!
  \file rd_efzc.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_efzc(): Reads in the frozen-core energy.
**
**   takes no arguments.
**
**   returns: double efzc  the frozen-core energy.
*/


double file30_rd_efzc(void)
{
  int  natom;
  PSI_FPTR junk;
  PSI_FPTR ptr;
  double efzc;

  natom = file30_rd_natom();
  ptr = (PSI_FPTR) 
           ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+4*sizeof(double));

  wreadw(info30_.filenum, (char *) &efzc, (int) sizeof(double),
	 ptr, &junk);

  return efzc;
}
