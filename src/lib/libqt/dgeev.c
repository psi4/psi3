
#if FCLINK==1
#define F_DGEEV dgeev_
#elif FCLINK==2
#define F_DGEEV dgeev
#else
#define F_DGEEV DGEEV
#endif


int C_DGEEV(int n, double **a, int lda,
  double *wr, double *wi, double **vl, int ldvl, double **vr,
  int ldvr, double *work, int lwork, int info)
{
  char jobvl, jobvr;
  jobvl = 'V';
  jobvr = 'V';
  F_DGEEV(&jobvl, &jobvr, &n, &(a[0][0]), &lda, &(wr[0]), &(wi[0]),
       &(vl[0][0]), &ldvl, &(vr[0][0]), &ldvr, &(work[0]), &lwork, &info);

  return info;
}

