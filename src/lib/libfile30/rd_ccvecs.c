#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*
** file30_rd_ccvecs():	Reads in a matrix, rows of which are ALPHA (ccvecs[0])
**			and BETA (ccvecs[1]) matrices of coupling coefficients 
**			for open shells stored in lower triangular form. 
**			Coupling coefficients are defined NOT as in 
**			C.C.J.Roothaan Rev. Mod. Phys. 32, 179 (1960) as it's 
**			stated in the manual pages for CSCF, but according to 
**			Pitzer (...) and are **different** from those in 
**			Yamaguchi, Osamura, Goddard, and Schaefer's book 
**			"Analytic Derivative Methods in Ab Initio Molecular 
**			Electronic Structure Theory".
**
**			The relationship between Pitzer's and Yamaguchi's 
**			conventions are follows :
**			ALPHA = 1-2*a , BETA = 1+4*b , where a and b are 
**			alpha's and beta's for open shells defined on pp. 69-70 
**			of Dr. Yamaguchi's book.
**
**   takes no arguments.
**
**   returns: 
**	double **ccvecs	a matrix 2 by abs(IOPEN) rows of which are coupling 
**			coefficient matrices for open-shells in packed form.
*/


double **file30_rd_ccvecs(void)
{
  int i, j;
  int nsymhf, nmo, mxcoef, ccvec_length;
  double **ccvecs;
  PSI_FPTR junk;
  PSI_FPTR alpha_ptr, beta_ptr;
  int tmp;
  int *scf_ptrs;
  
  nsymhf = file30_rd_nsymhf();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  ccvec_length = abs(file30_rd_iopen());
  
/* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  scf_ptrs = file30_rd_scf_ptrs();
  alpha_ptr = (PSI_FPTR) (scf_ptrs[8] -1)*sizeof(int);
  beta_ptr = (PSI_FPTR) (scf_ptrs[9] -1)*sizeof(int);
  free(scf_ptrs);
  
  if (ccvec_length > 0) {
    ccvecs = init_matrix(2,ccvec_length);

  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0]+60-1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr,
     &junk);
  alpha_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  alpha_ptr += (mxcoef+nmo)*sizeof(double) + (PSI_FPTR) (4*nsymhf*sizeof(int));*/
 
 /* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */

    wreadw(info30_.filenum, (char *) ccvecs[0], ccvec_length*sizeof(double),
	   alpha_ptr, &alpha_ptr);
    wreadw(info30_.filenum, (char *) ccvecs[1], ccvec_length*sizeof(double),
	   beta_ptr, &beta_ptr);

    return ccvecs;
  }
  else return NULL;
}
