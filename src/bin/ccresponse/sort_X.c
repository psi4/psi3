#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sort_X(char *cart, int irrep, double omega)
{
  dpdbuf4 X;
  char lbl[32];

  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart, omega);
  dpd_buf4_init(&X, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%1s_IAjb (%5.3f)", cart, omega);
  dpd_buf4_sort(&X, CC_LR, prqs, 10, 10, lbl);
  sprintf(lbl, "X_%1s_IbjA (%5.3f)", cart, omega);
  dpd_buf4_sort(&X, CC_LR, psqr, 10, 10, lbl);
  sprintf(lbl, "X_%1s_(2IjAb-IjbA) (%5.3f)", cart, omega);
  dpd_buf4_scmcopy(&X, CC_LR, lbl, 2);
  dpd_buf4_sort_axpy(&X, CC_LR, pqsr, 0, 5, lbl, -1);
  dpd_buf4_close(&X);

  sprintf(lbl, "X_%1s_IAjb (%5.3f)", cart, omega);
  dpd_buf4_init(&X, CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%1s_(2IAjb-IbjA) (%5.3f)", cart, omega);
  dpd_buf4_scmcopy(&X, CC_LR, lbl, 2);
  dpd_buf4_sort_axpy(&X, CC_LR, psrq, 10, 10, lbl, -1);
  dpd_buf4_close(&X);
}
