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

extern MOInfo_t MOInfo;
extern FILE *outfile;

extern void done(const char *);

double eval_rohf_derwfn_overlap()
{
  int nalpha = MOInfo.ndocc + MOInfo.nsocc;
  int nbeta = MOInfo.ndocc;
  double **CSC = eval_S_alpha();
  //    fprintf(outfile,"  -Cp*Spm*Cm :\n");
  //    print_mat(CSC,num_mo,num_mo,outfile);

  // Extract the alpha block
  double **CSC_alpha = block_matrix(nalpha,nalpha);
  for(int i=0;i<nalpha;i++)
    for(int j=0;j<nalpha;j++)
      CSC_alpha[i][j] = CSC[i][j];
  // Extract the beta block
  double **CSC_beta = block_matrix(nbeta,nbeta);
  for(int i=0;i<nbeta;i++)
    for(int j=0;j<nbeta;j++)
      CSC_beta[i][j] = CSC[i][j];
  free_block(CSC);

  // Compute the overlap of alpha part
  int *tmpintvec = new int[nalpha];
  C_DGETRF(nalpha,nalpha,&(CSC_alpha[0][0]),nalpha,tmpintvec);
  delete[] tmpintvec;
  double deter_a = 1.0;
  for(int i=0;i<nalpha;i++)
    deter_a *= CSC_alpha[i][i];
  deter_a = fabs(deter_a);
  free_block(CSC_alpha);

  // Compute the overlap of beta part
  tmpintvec = new int[nbeta];
  C_DGETRF(nbeta,nbeta,&(CSC_beta[0][0]),nbeta,tmpintvec);
  delete[] tmpintvec;
  double deter_b = 1.0;
  for(int i=0;i<nbeta;i++)
    deter_b *= CSC_beta[i][i];
  deter_b = fabs(deter_b);
  free_block(CSC_beta);

  return deter_a*deter_b;
}

