#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_beta_evals():  Writes the alpha SCF eigenvalues: the SCF orbital energies.
**
**  arguments: double *evals   an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
**
**  returns nothing
*/

void file30_wt_beta_evals(double *evals)
{
  PSI_FPTR scf_ptr, evalsptr, mo_coeff_ptr, junk;
  int *scf_ptrs;
  int tmp;
  
/* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */

  scf_ptrs = file30_rd_scf_ptrs();
  
  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 -1)*sizeof(int);*/

  /*wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);*/

  /*mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
    evalsptr = mo_coeff_ptr + (PSI_FPTR) file30_rd_mxcoef()*sizeof(double);*/

/* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
 
  evalsptr = (PSI_FPTR) (scf_ptrs[3] - 1)*sizeof(int);

  wwritw(info30_.filenum, (char *) evals, file30_rd_nmo()*sizeof(double),
	 evalsptr, &junk);

}
