#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double HXY(char *cart, int irrep, double omega)
{
  double polar;
  dpdfile2 X1, Y1, z;
  dpdbuf4 I;
  char lbl[32];

  sprintf(lbl, "Z_%1s_IA", cart);
  dpd_file2_init(&z, CC_TMP0, irrep, 0, 1, lbl);

  sprintf(lbl, "X_%1s_IA (-%5.3f)", cart, omega);
  dpd_file2_init(&Y1, CC_OEI, irrep, 0, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&Y1, &I, &z, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_file2_close(&Y1);

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  polar = 2.0 * dpd_file2_dot(&X1, &z);
  dpd_file2_close(&X1);

  dpd_file2_close(&z);

  return polar;
}
