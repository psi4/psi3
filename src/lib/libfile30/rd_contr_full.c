#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_contr_full(): Reads in the normalized contraction coefficients.
**
**  takes no arguments.
**
**  returns: double **contr
*/


double **file30_rd_contr_full(void)
{
  double **contr;
  double *temp_contr;
  int nprim, i, j, k, ij = 0;
  PSI_FPTR junk;
  PSI_FPTR contr_ptr;

  nprim = file30_rd_nprim();
  contr_ptr = (PSI_FPTR) (info30_.mpoint[5] - 1)*sizeof(int);

  temp_contr = init_array(MAXANGMOM*nprim);
  contr = block_matrix(nprim,MAXANGMOM);

  wreadw(info30_.filenum, (char *) temp_contr, (int) MAXANGMOM*nprim*sizeof(double),
	 contr_ptr, &junk);

  for(i=0;i<MAXANGMOM;i++) 
    for(k=0;k<nprim;k++) {
      contr[k][i] = temp_contr[ij];
      ij++;
    }

  free(temp_contr);

  return contr;
}
