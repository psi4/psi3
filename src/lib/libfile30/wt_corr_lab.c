/*!
  \file wt_corr_lab.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_corr_lab(): Writes a label into file30 which should describe the
**     wavefunction used to get the correlated energy to be stored in
**     file30 (see wt_ecor() description).
**
**   arguments:
** \param char *corr_lab  which should be a string like "CISD", or
**        "MCSCF" or some such wavefunction designation. Up to 7 character
**        + '\0'
**
**   returns: nothing.
**
**   N.B. The placement of the correlated energy in file30 is currently
**        under discussion, and this function may disappear with little or
**        no notice.
*/

void file30_wt_corr_lab(char *corrlab)
{
  int natom, tempi;
  PSI_FPTR junk;
  PSI_FPTR corrlab_ptr;

  natom = file30_rd_natom();
  tempi = info30_.mcalcs[0] + 60 + 20 + natom*6 - 1;
  corrlab_ptr = (PSI_FPTR) ((tempi)*sizeof(int) + 4*sizeof(double));

  wwritw(info30_.filenum, (char *) corrlab, sizeof(char)*8,
	 corrlab_ptr, &junk);
}
