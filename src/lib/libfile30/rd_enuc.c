/*!
  \file rd_enuc.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_enuc(): Reads in the nuclear repulsion energy
**
**   takes no arguments.
**
**   returns: double enuc  the nuclear repulsion energy.
*/

double file30_rd_enuc(void)
{
  int natom;
  PSI_FPTR junk;
  PSI_FPTR enuc_ptr;
  double enuc;

  natom = file30_rd_natom();
  enuc_ptr = (PSI_FPTR) (info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int);

  wreadw(info30_.filenum, (char *) &enuc, (int) sizeof(double),
	 enuc_ptr, &junk);

  return enuc;
}
