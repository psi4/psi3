#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double LCX(char *cart, int irrep, double omega)
{
  double polar=0.0;
  dpdfile2 X1, mu1, z1, l1;
  dpdbuf4 X2, mu2, z2, l2;
  char lbl[32];

  sprintf(lbl, "Mu_%1s_IA", cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  polar += 2.0 * dpd_file2_dot(&mu1, &X1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "z_%1s_IA", cart);
  dpd_file2_init(&z1, CC_TMP0, 0, 0, 1, lbl);

  sprintf(lbl, "MuBAR_%1s_MI", cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 0, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_contract222(&mu1, &X1, &z1, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "MuBAR_%1s_AE", cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 1, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_contract222(&X1, &mu1, &z1, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&mu1);

  sprintf(lbl, "MuBAR_%1s_ME", cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%1s_(2IjAb-IjbA) (%5.3f)", cart, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&mu1, &X2, &z1, 0, 0, 1, 1);
  dpd_buf4_close(&X2);
  dpd_file2_close(&mu1);

  dpd_file2_init(&l1, CC_OEI, 0, 0, 1, "LIA");
  polar += 2.0 * dpd_file2_dot(&z1, &l1);
  dpd_file2_close(&l1);

  dpd_file2_close(&z1);

  return polar;
}
