/*!
  \file rd_corr_lab.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_corr_lab(): Reads in a label from file30 which should describe the
**     wavefunction used to get the correlated energy which is stored in
**     file30 (see rd_ecorr() description).
**
**   takes no arguments.
**
**   returns a char *corr_lab  which should be a string like "CISD", or 
**        "MCSCF" or some such wavefunction designation. Up to 7 characters
**         + end of string character at the end
**   
**   N.B. The placement of the correlated energy in file30 is currently
**        under discussion, and this function may disappear with little or
**        no notice.
*/

char *file30_rd_corr_lab(void)
{
  int natom, tempi;
  PSI_FPTR junk;
  PSI_FPTR corrlab_ptr;
  char *corrlab;

  natom = file30_rd_natom();
  tempi = info30_.mcalcs[0] + 60 + 20 + natom*6 - 1;
  corrlab_ptr = (PSI_FPTR) ((tempi)*sizeof(int) + 4*sizeof(double));
 
  corrlab = (char *)malloc(sizeof(char)*8);

  wreadw(info30_.filenum, (char *) corrlab, 8*sizeof(char),
	 corrlab_ptr, &junk);
  corrlab[7] = '\0';

  return corrlab;
}
