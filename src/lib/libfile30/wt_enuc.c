/*!
  \file wt_enuc.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_wt_enuc(): Writes the nuclear repulsion energy
**
**   takes:  
** \param double enuc  the nuclear repulsion energy.
**
**   returns: none.
*/

void file30_wt_enuc(double enuc)
{
  int natom;
  PSI_FPTR junk;
  PSI_FPTR enuc_ptr;

  natom = file30_rd_natom();
  enuc_ptr = (PSI_FPTR) (info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int);

  wwritw(info30_.filenum, (char *) &enuc, (int) sizeof(double),
	 enuc_ptr, &junk);

  return;
}
