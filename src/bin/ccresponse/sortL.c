#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* sortL(): Sorts the MO-basis magnetic-dipole integrals into CC ordering
** for use in building the similarity-transformed integrals and certain
** components of the total linear response function.  
**
** NB: We must be careful when using the L integrals, because they are
** antisymmetric, i.e. LAI = -LIA.
**
** TDC, 6/03
*/

void sortL(void)
{
  int p, q, Gp, Gq, P, Q;
  int irrep;
  dpdfile2 l;
  double **L;

  /***** X-Component *****/
  irrep = moinfo.irrep_Rx;
  L = moinfo.LX;

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_X_IJ");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q]+moinfo.occ_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_X_AB");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p]+moinfo.vir_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  /***** Y-Component *****/
  irrep = moinfo.irrep_Ry;
  L = moinfo.LY;

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Y_IJ");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q]+moinfo.occ_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Y_AB");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p]+moinfo.vir_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  /***** Z-Component *****/
  irrep = moinfo.irrep_Rz;
  L = moinfo.LZ;

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Z_IJ");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.occpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_occ[q]+moinfo.occ_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Z_AB");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p]+moinfo.vir_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_file2_mat_init(&l);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p]+moinfo.occ_off[Gp]]-moinfo.nfzc;
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
	Q = moinfo.qt2pitzer[moinfo.qt_vir[q]+moinfo.vir_off[Gq]]-moinfo.nfzc;
	l.matrix[Gp][p][q] = L[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&l);
  dpd_file2_mat_close(&l);
  dpd_file2_close(&l);

}
