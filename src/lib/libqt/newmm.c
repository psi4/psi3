#include <stdio.h>

#if FCLINK==1
#define F_DGEMMC dgemmc_
#elif FCLINK==2
#define F_DGEMMC dgemmc
#else
#define F_DGEMMC DGEMMC
#endif

/*
** NEWMM
** 
** This is just a wrapper for the BLAS call DGEMM.  The only confusing
** thing is that C and Fortran have opposite stride conventions, so
** everything is transposed when the actual DGEMM call is made.
**
** This performs C = alpha*(opT)A * (opT)B + beta*C
*/

void newmm(double **A, int transa, double **B, int transb, double **C,
	   int num_rows, int num_links, int num_cols, 
           double alpha, double beta)
{
  int nra,nrb,nca,ncb;
  
  if(!transa) { nra = num_rows; nca = num_links; }
  else if(transa) {nra = num_links; nca = num_rows; }
  if(!transb) { nrb = num_links; ncb = num_cols; }
  else if(transb) { nrb = num_cols; ncb = num_links; }

  /* Error checking */
  if(!num_rows || !num_cols || !num_links) return;
  
  F_DGEMMC(&transa,&transb,&num_rows,&num_cols,&num_links,&alpha,
	   &(A[0][0]),&nra,&nca,&(B[0][0]),&nrb,&ncb,&beta,&(C[0][0]));
}
