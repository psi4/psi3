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

double eval_rhf_derwfn_overlap()
{
  int ndocc = MOInfo.ndocc;
  double **CSC = eval_S_alpha();
  //    fprintf(outfile,"  -Cp*Spm*Cm :\n");
  //    print_mat(CSC,num_mo,num_mo,outfile);

  // Extract the occupied block
  double **CSC_occ = block_matrix(ndocc,ndocc);
  for(int i=0;i<ndocc;i++)
    for(int j=0;j<ndocc;j++)
      CSC_occ[i][j] = CSC[i][j];
  free_block(CSC);

  // Compute the determinant
  int *tmpintvec = new int[ndocc];
  C_DGETRF(ndocc,ndocc,&(CSC_occ[0][0]),ndocc,tmpintvec);
  delete[] tmpintvec;
  double deter1 = 1.0;
  for(int i=0;i<ndocc;i++)
    deter1 *= CSC_occ[i][i];
  deter1 = fabs(deter1);
  //    fprintf(outfile,"  -Determinant for disp %d is %25.15lf\n\n",(disp-1)/2 + 1, deter1);

  free_block(CSC_occ);
  return deter1*deter1;
}

