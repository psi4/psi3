#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void denom1(char *cart, int irrep, double omega, dpdfile2 *X1);
void denom2(char *cart, int irrep, double omega, dpdbuf4 *X2);

void init_X(char *cart, int irrep, double omega)
{
  char lbl[32];
  dpdfile2 mu1, X1, FAE, FMI;
  dpdbuf4 X2, mu2;

  sprintf(lbl, "MuBAR_%1s_IA", cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_copy(&mu1, CC_OEI, lbl);
  dpd_file2_close(&mu1);

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  denom1(cart, irrep, omega, &X1);
  dpd_file2_close(&X1);

  sprintf(lbl, "MuBAR_%1s_IjAb", cart);
  dpd_buf4_init(&mu2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart, omega);
  dpd_buf4_copy(&mu2, CC_LR, lbl);
  dpd_buf4_close(&mu2);

  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  denom2(cart, irrep, omega, &X2);
  dpd_buf4_close(&X2);
}
