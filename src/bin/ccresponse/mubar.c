#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void mubar(void)
{
  int irrep;
  dpdfile2 mubar1, mu, t1, z;
  dpdbuf4 t2, mubar2;

  /***** X-component *****/

  irrep = moinfo.irrep_x;

  /** MuBAR_ME **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_X_ME");
  dpd_file2_close(&mu);

  /** MuBAR_MI **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_X_IJ");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_X_MI");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 0, "MuBAR_X_MI");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_AE **/
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_X_AB");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_X_AE");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 1, 1, "MuBAR_X_AE");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_IA **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_X_IA");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 1, "MuBAR_X_IA");
  if(!strcmp(params.gauge,"VELOCITY"))
    dpd_file2_scm(&mubar1, -1); /* the p integrals are antisymmetric! */

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_X_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_X_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &mu, &mubar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_X_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_contract222(&mu, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&mu);
  dpd_contract222(&z, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&mubar1);

  /** MuBAR_MbIj **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 10, 0, 10, 0, 0, "MuBAR_X_MbIj");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 5, 11, 5, 11, 0, "MuBAR_X_AbEi");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&mu, &t2, &mubar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_IjAb **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 0, 5, 0, 5, 0, "MuBAR_X_IjAb");

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_X_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 3, 1, 0, 1, 0);
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_X_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&mu, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_X_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_X_IA");
  dpd_contract222(&t1, &mu, &z, 1, 1, -1, 0);
  dpd_file2_close(&mu);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_X_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&mubar2);

  /***** Y-component *****/

  irrep = moinfo.irrep_y;

  /** MuBAR_ME **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Y_ME");
  dpd_file2_close(&mu);

  /** MuBAR_MI **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Y_IJ");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Y_MI");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 0, "MuBAR_Y_MI");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_AE **/
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Y_AB");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Y_AE");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 1, 1, "MuBAR_Y_AE");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_IA **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Y_IA");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 1, "MuBAR_Y_IA");
  if(!strcmp(params.gauge,"VELOCITY"))
    dpd_file2_scm(&mubar1, -1); /* the p integrals are antisymmetric! */

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Y_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Y_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &mu, &mubar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Y_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_contract222(&mu, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&mu);
  dpd_contract222(&z, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&mubar1);

  /** MuBAR_MbIj **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 10, 0, 10, 0, 0, "MuBAR_Y_MbIj");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 5, 11, 5, 11, 0, "MuBAR_Y_AbEi");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&mu, &t2, &mubar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_IjAb **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 0, 5, 0, 5, 0, "MuBAR_Y_IjAb");

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Y_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 3, 1, 0, 1, 0);
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Y_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&mu, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_Y_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Y_IA");
  dpd_contract222(&t1, &mu, &z, 1, 1, -1, 0);
  dpd_file2_close(&mu);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Y_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&mubar2);

  /***** Z-component *****/

  irrep = moinfo.irrep_z;

  /** MuBAR_ME **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Z_ME");
  dpd_file2_close(&mu);

  /** MuBAR_MI **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Z_IJ");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Z_MI");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 0, "MuBAR_Z_MI");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_AE **/
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Z_AB");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Z_AE");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 1, 1, "MuBAR_Z_AE");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);
  dpd_file2_close(&mubar1);

  /** MuBAR_IA **/
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_file2_copy(&mu, CC_OEI, "MuBAR_Z_IA");
  dpd_file2_close(&mu);

  dpd_file2_init(&mubar1, CC_OEI, irrep, 0, 1, "MuBAR_Z_IA");
  if(!strcmp(params.gauge,"VELOCITY"))
    dpd_file2_scm(&mubar1, -1); /* the p integrals are antisymmetric! */

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Z_AB");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&t1, &mu, &mubar1, 0, 0, 1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Z_IJ");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&mu, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&t2, &mu, &mubar1, 0, 0, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Z_MI");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_contract222(&mu, &t1, &z, 0, 0, 1, 0);
  dpd_file2_close(&mu);
  dpd_contract222(&z, &t1, &mubar1, 1, 1, -1, 1);
  dpd_file2_close(&t1);
  dpd_file2_close(&z);

  dpd_file2_close(&mubar1);

  /** MuBAR_MbIj **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 10, 0, 10, 0, 0, "MuBAR_Z_MbIj");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 0, 1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_AbEi -- stored (Ab,Ei) **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 5, 11, 5, 11, 0, "MuBAR_Z_AbEi");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");
  dpd_contract244(&mu, &t2, &mubar2, 0, 2, 1, -1, 0);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);
  dpd_buf4_close(&mubar2);

  /** MuBAR_IjAb **/
  dpd_buf4_init(&mubar2, CC_LR, irrep, 0, 5, 0, 5, 0, "MuBAR_Z_IjAb");

  dpd_file2_init(&mu, CC_OEI, irrep, 1, 1, "Mu_Z_AB");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 3, 1, 0, 1, 0);
  dpd_contract244(&mu, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&mu, CC_OEI, irrep, 0, 0, "Mu_Z_IJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &mu, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&mu, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&mu);

  dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, "z_Z_AE");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&mu, CC_OEI, irrep, 0, 1, "Mu_Z_IA");
  dpd_contract222(&t1, &mu, &z, 1, 1, -1, 0);
  dpd_file2_close(&mu);
  dpd_file2_close(&t1);
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 3, 1, 0, 1, 1);
  dpd_contract244(&z, &t2, &mubar2, 1, 2, 1, 1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, "z_Z_MI"); /* generated above */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&t2, &z, &mubar2, 1, 0, 1, -1, 1);
  dpd_contract244(&z, &t2, &mubar2, 0, 0, 0, -1, 1);
  dpd_buf4_close(&t2);
  dpd_file2_close(&z);

  dpd_buf4_close(&mubar2);
}
