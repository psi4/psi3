/*!
  \file rd_alpha_lagr.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_alpha_lagr():	Reads in an alpha lagrangian matrix in MO basis
**
**   takes no arguments.
**
**   returns: 
**	double **lagr	a matrix nmo by nmo.
*/


double **file30_rd_alpha_lagr(void)
{
  int i, j;
  int nsymhf, nmo, mxcoef, iopen;
  double *lagr_tri, **lagr_sq;
  PSI_FPTR lagr_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;

  nsymhf = file30_rd_nsymhf();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  iopen = file30_rd_iopen();

  /* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  
  scf_ptrs = file30_rd_scf_ptrs();
  
  lagr_tri = init_array(nmo*(nmo+1)/2);
  
  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0]+60-1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr,
     &scf_ptr);
  lagr_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  lagr_ptr += (mxcoef+nmo)*sizeof(double) +
  (PSI_FPTR) (nsymhf*4*sizeof(char) + 2*nsymhf*sizeof(int));
  
  if (iopen != 0)
    lagr_ptr += (PSI_FPTR) (nsymhf*sizeof(int) + 2*abs(iopen)*sizeof(double));*/

/* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  lagr_ptr = (PSI_FPTR) (scf_ptrs[10]-1)*sizeof(int);
  
  wreadw(info30_.filenum, (char *) lagr_tri, nmo*(nmo+1)*sizeof(double)/2,
               lagr_ptr, &lagr_ptr);

  lagr_sq = block_matrix(nmo,nmo);
  tri_to_sq(lagr_tri,lagr_sq,nmo);
  free(lagr_tri);

  free(scf_ptrs);

  return lagr_sq;
  
}
