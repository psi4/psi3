/*!
  \file wt_iopen.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_iopen():  Writes out the total number of open shells.
**
**   arguments:
** \param int iopen: the total number of singly occupied orbitals in the SCF
**     reference wavefunction.
**
**   returns: nothing
**
*/


void file30_wt_iopen(int iopen)
{
  PSI_FPTR junk;
  
  info30_.mconst[42] = iopen;

  /* Write out the complete list of constants to file30 */
  wwritw(info30_.filenum, (char *) info30_.mconst,
	 sizeof(unsigned int)*MCONST, sizeof(unsigned int)*100,
	 &junk);
}
