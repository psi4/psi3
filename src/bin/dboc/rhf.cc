#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
}
#include "params.h"
#include "moinfo.h"
#include "mo_overlap.h"
#include "float.h"
#include "linalg.h"

extern MOInfo_t MOInfo;
extern FILE *outfile;
extern Params_t Params;

extern void done(const char *);

double eval_rhf_derwfn_overlap()
{
  int ndocc = MOInfo.ndocc;
  FLOAT **CSC = eval_S_alpha();

  chkpt_init(PSIO_OPEN_OLD);
  int* clsdpi = chkpt_rd_clsdpi();
  int* orbspi = chkpt_rd_orbspi();
  int nirreps = chkpt_rd_nirreps();
  chkpt_close();

  // Extract the occupied block
  FLOAT **CSC_occ = create_matrix(ndocc,ndocc);
  int mo_offset1 = 0;
  int occ_offset1 = 0;
  for(int irrep1=0; irrep1<nirreps; irrep1++) {

    int nocc1 = clsdpi[irrep1];

    int mo_offset2 = 0;
    int occ_offset2 = 0;
    for(int irrep2=0; irrep2<nirreps; irrep2++) {

      int nocc2 = clsdpi[irrep2];

      for(int i=0;i<nocc1;i++)
	for(int j=0;j<nocc2;j++)
	  CSC_occ[i+occ_offset1][j+occ_offset2] = CSC[i+mo_offset1][j+mo_offset2];

      occ_offset2 += nocc2;
      mo_offset2 += orbspi[irrep2];
    }

    occ_offset1 += nocc1;
    mo_offset1 += orbspi[irrep1];
  }
  delete_matrix(CSC);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile,"  +/- overlap in the basis of doubly-occupied MOs:\n");
    print_mat(CSC_occ, ndocc, ndocc, outfile);
  }

  // Compute the determinant
  int *tmpintvec = new int[ndocc];
  FLOAT sign;
  lu_decom(CSC_occ, ndocc, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter1 = 1.0;
  for(int i=0;i<ndocc;i++)
    deter1 *= CSC_occ[i][i];
  deter1 = FABS(deter1);

  delete_matrix(CSC_occ);
  return (double)deter1*deter1;
}

