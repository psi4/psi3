#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double LCX(char *cart_c, int irrep_c, char *cart_x, int irrep_x, double omega)
{
  double polar=0.0;
  dpdfile2 X1, mu1, z1, l1;
  dpdbuf4 X2, mu2, z2, l2, Z;
  char lbl[32];

  /*** Mu * X1 ***/

  sprintf(lbl, "Mu_%1s_IA", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  polar += 2.0 * dpd_file2_dot(&mu1, &X1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  /*** L1 * MuBAR * X1 + L1 * MuBAR * X2 ***/

  dpd_file2_init(&z1, CC_TMP0, 0, 0, 1, "z_IA");

  sprintf(lbl, "MuBAR_%1s_MI", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 0, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&mu1, &X1, &z1, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "MuBAR_%1s_AE", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 1, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &mu1, &z1, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "MuBAR_%1s_ME", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%1s_(2IjAb-IjbA) (%5.3f)", cart_x, omega);
  dpd_buf4_init(&X2, CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&mu1, &X2, &z1, 0, 0, 1, 1);
  dpd_buf4_close(&X2);
  dpd_file2_close(&mu1);

  dpd_file2_init(&l1, CC_OEI, 0, 0, 1, "LIA");
  polar += 2.0 * dpd_file2_dot(&z1, &l1);
  dpd_file2_close(&l1);

  dpd_file2_close(&z1);

  /*** L2 * MuBAR * X1 + L2 * MuBAR * X2 ***/

  dpd_buf4_init(&z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&z2, 0);

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "MuBAR_%1s_EiAb", cart_c);
  dpd_buf4_init(&mu2, CC_LR, irrep_c, 11, 5, 11, 5, 0, lbl);
  dpd_contract244(&X1, &mu2, &Z, 1, 0, 0, 1, 0);
  dpd_buf4_close(&mu2);
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "MuBAR_%1s_MbIj", cart_c);
  dpd_buf4_init(&mu2, CC_LR, irrep_c, 10, 0, 10, 0, 0, lbl);
  dpd_contract244(&X1, &mu2, &Z, 0, 0, 1, 1, 0);
  dpd_buf4_close(&mu2);
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_close(&Z);

  dpd_file2_close(&X1);

  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_x, omega);
  dpd_buf4_init(&X2, CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "MuBAR_%1s_AE", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 1, 1, lbl);
  dpd_contract424(&X2, &mu1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&mu1);
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, 1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "MuBAR_%1s_MI", cart_c);
  dpd_file2_init(&mu1, CC_OEI, irrep_c, 0, 0, lbl);
  dpd_contract244(&mu1, &X2, &Z, 0, 0, 0, 1, 0);
  dpd_file2_close(&mu1);
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_sort(&Z, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z, &z2, -1);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&X2);

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += dpd_buf4_dot(&l2, &z2);
  dpd_buf4_close(&l2);

  dpd_buf4_close(&z2);

  return polar;
}
