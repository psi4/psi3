#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_ecorr():  Writes the correlated energy to file30.
**
**    arguments: double ecorr   the correlated energy.  To write some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see wt_corr_lab().
**
**   N.B. The placement of the correlated energy in file30 is currently
**        under discussion, and this function may disappear with little or
**        no notice.
*/

void file30_wt_ecorr(double ecorr)
{
  int natom;
  PSI_FPTR junk;
  PSI_FPTR ecorr_ptr;

  natom = file30_rd_natom();
  ecorr_ptr = (PSI_FPTR) 
        ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+3*sizeof(double));

  wwritw(info30_.filenum, (char *) &ecorr, (int) sizeof(double),
	 ecorr_ptr, &junk);

}
