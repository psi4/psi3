#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_ecorr():  Reads in the correlated energy stored in file30.
**
**    takes no arguments.
**
**    returns: double e_corr   the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
**   
**   N.B. The placement of the correlated energy in file30 is currently
**        under discussion, and this function may disappear with little or
**        no notice.
*/

double file30_rd_ecorr(void)
{
  int natom;
  PSI_FPTR junk;
  PSI_FPTR ecorr_ptr;
  double ecorr;

  natom = file30_rd_natom();
  ecorr_ptr = (PSI_FPTR) 
        ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+6*sizeof(double));

  wreadw(info30_.filenum, (char *) &ecorr, (int) sizeof(double),
	 ecorr_ptr, &junk);

  return ecorr;
}
