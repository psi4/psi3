/*!
  \file rd_escf.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_escf(): Reads in the scf energy.
**
**   takes no arguments.
**
**   returns: double escf  the scf energy.
*/


double file30_rd_escf(void)
{
  int  natom;
  PSI_FPTR junk;
  PSI_FPTR escf_ptr;
  double escf;

  natom = file30_rd_natom();
  escf_ptr = (PSI_FPTR) 
           ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+sizeof(double));

  wreadw(info30_.filenum, (char *) &escf, (int) sizeof(double),
	 escf_ptr, &junk);

  return escf;
}
