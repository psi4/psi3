#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
}
#include "moinfo.h"
#include "mo_overlap.h"
#include "float.h"
#include "linalg.h"

extern MOInfo_t MOInfo;
extern FILE *outfile;

extern void done(const char *);

double eval_uhf_derwfn_overlap()
{
  int nalpha = MOInfo.nalpha;
  int nbeta = MOInfo.nbeta;
  FLOAT **CSC_a = eval_S_alpha();
  FLOAT **CSC_b = eval_S_beta();

  // Extract the alpha block
  FLOAT **CSC_alpha = create_matrix(nalpha,nalpha);
  for(int i=0;i<nalpha;i++)
    for(int j=0;j<nalpha;j++)
      CSC_alpha[i][j] = CSC_a[i][j];
  delete_matrix(CSC_a);
  // Extract the beta block
  FLOAT **CSC_beta = create_matrix(nbeta,nbeta);
  for(int i=0;i<nbeta;i++)
    for(int j=0;j<nbeta;j++)
      CSC_beta[i][j] = CSC_b[i][j];
  delete_matrix(CSC_b);

  // Compute the overlap of alpha part
  int *tmpintvec = new int[nalpha];
  //  C_DGETRF(nalpha,nalpha,&(CSC_alpha[0][0]),nalpha,tmpintvec);
  FLOAT sign;
  lu_decom(CSC_alpha, nalpha, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter_a = 1.0;
  for(int i=0;i<nalpha;i++)
    deter_a *= CSC_alpha[i][i];
  deter_a = FABS(deter_a);
  delete_matrix(CSC_alpha);

  // Compute the overlap of beta part
  tmpintvec = new int[nbeta];
  //  C_DGETRF(nbeta,nbeta,&(CSC_beta[0][0]),nbeta,tmpintvec);
  lu_decom(CSC_beta, nbeta, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter_b = 1.0;
  for(int i=0;i<nbeta;i++)
    deter_b *= CSC_beta[i][i];
  deter_b = FABS(deter_b);
  delete_matrix(CSC_beta);

  return (double)deter_a*deter_b;
}

