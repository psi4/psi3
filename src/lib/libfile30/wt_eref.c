/*!
  \file wt_eref.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_eref(): Writes the reference energy to file30.
**
**   arguments:
** \param double eref  the reference energy.
**
**   returns nothing.
*/

void file30_wt_eref(double eref)
{
  int  natom;
  PSI_FPTR junk;
  PSI_FPTR eref_ptr;

  natom = file30_rd_natom();
  eref_ptr = (PSI_FPTR) 
           ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+2*sizeof(double));

  wwritw(info30_.filenum, (char *) &eref, (int) sizeof(double),
	 eref_ptr, &junk);

  return;
}
