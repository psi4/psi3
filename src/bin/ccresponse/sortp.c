#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sortp(void)
{
  int p, q, Gp, Gq, P, Q;
  int irrep;
  dpdfile2 mu;
  double **MU;

  /***** X-Component *****/
  irrep = moinfo.irrep_x;
  MU = moinfo.MUX;

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_X_IJ");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q+moinfo.occ_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_X_AB");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p+moinfo.vir_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  /***** Y-Component *****/
  irrep = moinfo.irrep_y;
  MU = moinfo.MUY;

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Y_IJ");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q+moinfo.occ_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Y_AB");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p+moinfo.vir_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  /***** Z-Component *****/
  irrep = moinfo.irrep_z;
  MU = moinfo.MUZ;

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Z_IJ");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q+moinfo.occ_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Z_AB");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p+moinfo.vir_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_file2_mat_init(&mu);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	mu.matrix[Gp][p][q] = MU[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);
}
