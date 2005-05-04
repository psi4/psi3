#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* Lbar(): Constructs the components of the similarity-transform
** magnetic-dipole operator. 
**
** NB: We must be careful in contractions involving the bare
** magnetic-dipole integrals, because they are antisymmetric, i.e. LAI
** = - LIA.  This only affects the leading term of LBAR_IA (which is
** poorly named and should really be called LBAR_AI), where I copy the
** LIA integrals into place with a -1 sign.
**
** TDC, 6/03
*/

void Lbar(void)
{
  int irrep;
  dpdfile2 lbar1, l, t1, z;
  dpdbuf4 t2, lbar2;

  /***** X-component *****/

  irrep = moinfo.irrep_Rx;

  /** LBAR_ME **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_X_ME");
  dpd_file2_close(&l);

  /** LBAR_MI **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_X_IJ");
  dpd_file2_copy(&l, CC_OEI, "LBAR_X_MI");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 0, "LBAR_X_MI");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_AE **/
  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_X_AB");
  dpd_file2_copy(&l, CC_OEI, "LBAR_X_AE");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 1, 1, "LBAR_X_AE");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_IA **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_X_IA");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 1, "LBAR_X_IA");
  dpd_file2_scm(&lbar1, -1); /* the L integrals are antisymmetric! */

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_X_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_X_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &l, &lbar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_X_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_contract222(&l, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&l);
  dpd_contract222(&z, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&lbar1);

  /** LBAR_MbIj **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 10, 0, 10, 0, 0, "LBAR_X_MbIj");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 5, 11, 5, 11, 0, "LBAR_X_AbEi");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&l, &t2, &lbar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_IjAb **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 0, 5, 0, 5, 0, "LBAR_X_IjAb");

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_X_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 3, 1, 0, 1, 0);
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_X_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&l, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_X_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_X_IA");
  dpd_contract222(&t1, &l, &z, 1, 1, -1, 0);
  dpd_file2_close(&l);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_X_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&lbar2);

  /***** Y-component *****/

  irrep = moinfo.irrep_Ry;

  /** LBAR_ME **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Y_ME");
  dpd_file2_close(&l);

  /** LBAR_MI **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Y_IJ");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Y_MI");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 0, "LBAR_Y_MI");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_AE **/
  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Y_AB");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Y_AE");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 1, 1, "LBAR_Y_AE");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_IA **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Y_IA");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 1, "LBAR_Y_IA");
  dpd_file2_scm(&lbar1, -1); /* the L integrals are antisymmetric! */

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Y_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Y_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &l, &lbar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Y_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_contract222(&l, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&l);
  dpd_contract222(&z, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&lbar1);

  /** LBAR_MbIj **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 10, 0, 10, 0, 0, "LBAR_Y_MbIj");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 5, 11, 5, 11, 0, "LBAR_Y_AbEi");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&l, &t2, &lbar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_IjAb **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 0, 5, 0, 5, 0, "LBAR_Y_IjAb");

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Y_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 3, 1, 0, 1, 0);
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Y_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&l, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_Y_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Y_IA");
  dpd_contract222(&t1, &l, &z, 1, 1, -1, 0);
  dpd_file2_close(&l);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Y_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&lbar2);

  /***** Z-component *****/

  irrep = moinfo.irrep_Rz;

  /** LBAR_ME **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Z_ME");
  dpd_file2_close(&l);

  /** LBAR_MI **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Z_IJ");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Z_MI");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 0, "LBAR_Z_MI");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_AE **/
  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Z_AB");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Z_AE");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 1, 1, "LBAR_Z_AE");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);
  dpd_file2_close(&lbar1);

  /** LBAR_IA **/
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_file2_copy(&l, CC_OEI, "LBAR_Z_IA");
  dpd_file2_close(&l);

  dpd_file2_init(&lbar1, CC_OEI, irrep, 0, 1, "LBAR_Z_IA");
  dpd_file2_scm(&lbar1, -1); /* the L integrals are antisymmetric! */

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Z_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &l, &lbar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Z_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&l, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &l, &lbar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Z_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_contract222(&l, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&l);
  dpd_contract222(&z, &t1, &lbar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&lbar1);

  /** LBAR_MbIj **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 10, 0, 10, 0, 0, "LBAR_Z_MbIj");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 5, 11, 5, 11, 0, "LBAR_Z_AbEi");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&l, &t2, &lbar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);
  dpd_buf4_close(&lbar2);

  /** LBAR_IjAb **/
  dpd_buf4_init(&lbar2, CC_LR, irrep, 0, 5, 0, 5, 0, "LBAR_Z_IjAb");

  dpd_file2_init(&l, CC_OEI, irrep, 1, 1, "L_Z_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 3, 1, 0, 1, 0);
  dpd_contract244(&l, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&l, CC_OEI, irrep, 0, 0, "L_Z_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &l, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&l, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&l);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_Z_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&l, CC_OEI, irrep, 0, 1, "L_Z_IA");
  dpd_contract222(&t1, &l, &z, 1, 1, -1, 0);
  dpd_file2_close(&l);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &lbar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Z_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &lbar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &lbar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&lbar2);

}
