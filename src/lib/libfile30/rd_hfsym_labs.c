#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_hfsym_labs(): Read in the symmetry labels _only_ for those irreps 
**     which have basis functions.  
**
**   takes no arguments.
**
**   returns: char **hfsym_labs   an array of labels (strings) which denote
**      the irreps which have basis functions (in Cotton ordering).  For DZ or
**      STO water, for example, in C2v symmetry, this would be an array of 
**      three labels: "A1", "B1", and "B2".
*/

char **file30_rd_hfsym_labs(void)
{
  PSI_FPTR scf_ptr, lab_ptr, mo_coeff_ptr, junk;
  int i, nsymhf, mxcoef, nmo, tmp;
  char **hfsym_labs;
  int *scf_ptrs;
  
  nsymhf = file30_rd_nsymhf();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  scf_ptrs = file30_rd_scf_ptrs();
  
/*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 -1)*sizeof(int);
  
  wreadw(info30_.filenum, (char *) scf_points, sizeof(int)*20,scf_ptr, &junk);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);*/

  /*mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);*/
  /*lab_ptr = mo_coeff_ptr + (PSI_FPTR) mxcoef*sizeof(double)+(PSI_FPTR) nmo*sizeof(double);*/

  lab_ptr = (PSI_FPTR) (scf_ptrs[4] - 1)*sizeof(int);
  hfsym_labs = (char **)malloc(sizeof(char *)*nsymhf);
  for(i=0;i<nsymhf;i++)
    {
     hfsym_labs[i] = (char *)malloc(sizeof(int));
     wreadw(info30_.filenum, (char *) hfsym_labs[i], sizeof(int),
	  lab_ptr, &lab_ptr);
     hfsym_labs[i][3] = '\0';
    }
  free(scf_ptrs);

  return hfsym_labs;

}
