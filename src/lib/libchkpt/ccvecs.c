/*!
  \file ccvecs.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_ccvecs():	Reads in a matrix, rows of which are ALPHA (ccvecs[0])
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
**			ALPHA = 1-2*a , BETA = -1-4*b , where a and b are 
**			alpha's and beta's for open shells defined on pp. 69-70 
**			of Dr. Yamaguchi's book.
**
**   takes no arguments.
**
**   returns: 
**	double **ccvecs	a matrix 2 by abs(IOPEN) rows of which are coupling 
**        coefficient matrices for open-shells in packed form.
** \ingroup (CHKPT)
*/
double **chkpt_rd_ccvecs(void)
{
  int nmo, ccvec_length;
  double **ccvecs;
  char *keyword;
  
  nmo = chkpt_rd_nmo();
  ccvec_length = abs(chkpt_rd_iopen());
  
  if (ccvec_length > 0) {
    ccvecs = block_matrix(2,ccvec_length);

    keyword = chkpt_build_keyword("SCF coupling coefficients"); 

    psio_read_entry(PSIF_CHKPT, keyword, (char *) ccvecs[0],
        2*ccvec_length*sizeof(double));

    free(keyword);

    return ccvecs;
  }
  else return NULL;
}


/*!
** chkpt_wt_ccvecs():	Writes a matrix of coupling coefficients.  See the 
**                      comments chkpt_rd_ccvecs() above.
**
**    \param ccvecs =   a matrix 2 by abs(IOPEN) rows of which are coupling 
**			coefficient matrices for open-shells in packed form.
**
** returns: none
** \ingroup (CHKPT)
*/
void chkpt_wt_ccvecs(double **ccvecs)
{
  int nmo, ccvec_length;
  char *keyword;
  
  nmo = chkpt_rd_nmo();
  ccvec_length = abs(chkpt_rd_iopen());
  
  if (ccvec_length > 0) {
    keyword = chkpt_build_keyword("SCF coupling coefficients");

    psio_write_entry(PSIF_CHKPT, keyword, (char *) ccvecs[0],
       2*ccvec_length*sizeof(double));

    free(keyword);
  }
}

