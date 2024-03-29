/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* SORTI_UHF(): Place all the components of the UHF Lagrangian into a
** large matrix, I (moinfo.I), which we also symmetrize by computing
** Ipq = 1/2 (Ipq + Iqp).  This matrix is later written to disk in
** dump() for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_UHF(void)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off; 
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off; 
  int *aocc_sym, *avir_sym;
  int *bocc_sym, *bvir_sym;
  int *qt_aocc, *qt_avir;
  int *qt_bocc, *qt_bvir;
  double **O_a, **O_b;
  double chksum, value;
  dpdfile2 D;

  nmo = moinfo.nmo;
  nfzc = moinfo.nfzc;
  nfzv = moinfo.nfzv;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;
  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi; bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off; avir_off = moinfo.avir_off;
  bocc_off = moinfo.bocc_off; bvir_off = moinfo.bvir_off;
  aocc_sym = moinfo.aocc_sym; avir_sym = moinfo.avir_sym;
  bocc_sym = moinfo.bocc_sym; bvir_sym = moinfo.bvir_sym;

  qt_aocc = moinfo.qt_aocc; qt_avir = moinfo.qt_avir;
  qt_bocc = moinfo.qt_bocc; qt_bvir = moinfo.qt_bvir;

  O_a = block_matrix(nmo,nmo);
  O_b = block_matrix(nmo,nmo);

  /* Sort alpha components first */
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "I(I,J)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(j=0; j < aoccpi[h]; j++) {
	J = qt_aocc[aocc_off[h] + j];
	O_a[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "I'AB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < avirtpi[h]; a++) {
      A = qt_avir[avir_off[h] + a];
      for(b=0; b < avirtpi[h]; b++) {
	B = qt_avir[avir_off[h] + b];

	O_a[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "I(I,A)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(a=0; a < avirtpi[h]; a++) {
	A = qt_avir[avir_off[h] + a];

	O_a[A][I] += D.matrix[h][i][a];
	O_a[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Sort beta components */
  dpd_file2_init(&D, CC_OEI, 0, 2, 2, "I(i,j)");
  dpd_file2_mat_init(&D); 
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) { 
      I = qt_bocc[bocc_off[h] + i];
      for(j=0; j < boccpi[h]; j++) {
	J = qt_bocc[bocc_off[h] + j];
	O_b[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 3, 3, "I'ab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < bvirtpi[h]; a++) {
      A = qt_bvir[bvir_off[h] + a];
      for(b=0; b < bvirtpi[h]; b++) {
	B = qt_bvir[bvir_off[h] + b];

	O_b[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 2, 3, "I(i,a)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      I = qt_bocc[bocc_off[h] + i];
      for(a=0; a < bvirtpi[h]; a++) {
	A = qt_bvir[bvir_off[h] + a];

	O_b[A][I] += D.matrix[h][i][a];
	O_b[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Symmetrize the Lagrangians */
  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5*(O_a[p][q] + O_a[q][p]);
      O_a[p][q] = O_a[q][p] = value;

      value = 0.5*(O_b[p][q] + O_b[q][p]);
      O_b[p][q] = O_b[q][p] = value;
    }
  }

  /*
  fprintf(outfile, "\n\tAlpha MO Lag:\n");
  mat_print(O_a, nmo, nmo, outfile);
  fprintf(outfile, "\n\tBeta MO Lag:\n");
  mat_print(O_b, nmo, nmo, outfile);
  */

  /* Multiply the Lagrangian by -2.0 for the final energy derivative
     expression */
  for(p=0; p < nmo; p++) {
    for(q=0; q < nmo; q++) {
      O_a[p][q] *= -2.0;
      O_b[p][q] *= -2.0;
    }
  }

  moinfo.I_a = O_a;
  moinfo.I_b = O_b;
}

}} // namespace psi::ccdensity
