/*!
  \file wt_disp.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_wt_disp():  Writes out the current displacement number.
**
** This routine is called by ugeom to update the displacement number.
**
**   arguments:
** \param int disp: the current geometry displacement number.
**
**   returns: nothing
**
*/


void file30_wt_disp(int disp)
{
  PSI_FPTR junk;
  
  info30_.mconst[51] = disp;

  /* Write out the complete list of constants to file30 */
  wwritw(info30_.filenum, (char *) info30_.mconst,
	 sizeof(unsigned int)*MCONST, sizeof(unsigned int)*100,
	 &junk);
}
