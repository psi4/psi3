/*!
  \file rd_eref.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_eref(): Reads in the reference energy (if different from SCF).
**
**   takes no arguments.
**
**   returns: double eref  the reference energy.
*/


double file30_rd_eref(void)
{
  int  natom;
  PSI_FPTR junk;
  PSI_FPTR eref_ptr;
  double eref;

  natom = file30_rd_natom();
  eref_ptr = (PSI_FPTR) 
           ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+2*sizeof(double));

  wreadw(info30_.filenum, (char *) &eref, (int) sizeof(double),
	 eref_ptr, &junk);

  return eref;
}
