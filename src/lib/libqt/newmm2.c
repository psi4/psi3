#include <stdio.h>

#if FCLINK==1
#define F_DGEMMC dgemmc_
#elif FCLINK==2
#define F_DGEMMC dgemmc
#else
#define F_DGEMMC DGEMMC
#endif

void newmm2(double *A, int transa, double *B, int transb, double *C,
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
	   &(A[0]),&nra,&nca,&(B[0]),&nrb,&ncb,&beta,&(C[0]));
}
