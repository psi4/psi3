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

double eval_rhf_derwfn_overlap()
{
  int ndocc = MOInfo.ndocc;
  FLOAT **CSC = eval_S_alpha();
  //    fprintf(outfile,"  -Cp*Spm*Cm :\n");
  //    print_mat(CSC,num_mo,num_mo,outfile);

  // Extract the occupied block
  FLOAT **CSC_occ = create_matrix(ndocc,ndocc);
  for(int i=0;i<ndocc;i++)
    for(int j=0;j<ndocc;j++)
      CSC_occ[i][j] = CSC[i][j];
  delete_matrix(CSC);

  // Compute the determinant
  int *tmpintvec = new int[ndocc];
  //  C_DGETRF(ndocc,ndocc,&(CSC_occ[0][0]),ndocc,tmpintvec);
  FLOAT sign;
  lu_decom(CSC_occ, ndocc, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter1 = 1.0;
  for(int i=0;i<ndocc;i++)
    deter1 *= CSC_occ[i][i];
  deter1 = FABS(deter1);
#if LONG_DOUBLE
  fprintf(outfile,"  -Determinant for disp %d is %25.15Lf\n\n",1, deter1);
#else
  //  fprintf(outfile,"  -Determinant for disp %d is %25.15lf\n\n",(disp-1)/2 + 1, deter1);
#endif

  delete_matrix(CSC_occ);
  return deter1*deter1;
}

