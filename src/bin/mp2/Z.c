#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void rhf_Z(void);
void uhf_Z(void);

void Z(void)
{
  if(params.ref == 0) rhf_Z();
  else if(params.ref == 2) uhf_Z();
}

void rhf_Z(void)
{
  dpdfile2 L;
  dpdfile2 D;
  dpdbuf4 A;
  double **Z;
  int h, nirreps;
  int a, i, num_ai, count;
  int I, B;

  nirreps = mo.nirreps;

  dpd_file2_init(&L, CC_OEI, 0, 1, 0, "LAI");
  dpd_file2_mat_init(&L);
  dpd_file2_mat_rd(&L);
  num_ai = 0;
  for(h=0; h < nirreps; h++)
    num_ai += L.params->rowtot[h]*L.params->coltot[h];

  Z = block_matrix(1,num_ai);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < L.params->rowtot[h]; a++)
      for(i=0; i < L.params->coltot[h]; i++) 
	Z[0][count++] = -L.matrix[h][a][i];

  dpd_file2_mat_close(&L);
  dpd_file2_close(&L);

  dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_mat_irrep_init(&A, 0);
  dpd_buf4_mat_irrep_rd(&A, 0);

  pople(A.matrix[0], Z[0], num_ai, 1, 1e-12, outfile, 0);

  dpd_buf4_mat_irrep_close(&A, 0);
  dpd_buf4_close(&A);

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
  dpd_file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) 
	D.matrix[h][a][i] = Z[0][count++];

  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /*
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  for(h=0; h<mo.nirreps; h++) {
    for(a=0; a<mo.virpi[h]; a++) {
      B = mo.qt_vir[mo.vir_off[h] + a];
      for(i=0; i<mo.occpi[h]; i++) {
        I = mo.qt_occ[mo.occ_off[h] + i];
        mo.opdm[B][I] += D.matrix[h][a][i];
        mo.opdm[I][B] += D.matrix[h][a][i];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
  */

  free_block(Z);
}

void uhf_Z(void)
{

}
