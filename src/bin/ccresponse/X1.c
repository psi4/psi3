#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void denom1(dpdfile2 *X1, double omega);
void local_filter_T1(dpdfile2 *T1, double omega);

void X1_build(char *pert, char *cart, int irrep, double omega)
{
  dpdfile2 F, X1, X1new;
  dpdbuf4 W, X2;
  char lbl[32];

  sprintf(lbl, "%sBAR_%1s_IA", pert, cart);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "New X_%s_%1s_IA (%5.3f)", pert, cart, omega);
  dpd_file2_copy(&X1new, CC_OEI, lbl);
  dpd_file2_close(&X1new);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);

  /*** S-S ***/

  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert, cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);

  dpd_file2_axpy(&X1, &X1new, -omega, 0);

  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "FAE");
  dpd_contract222(&X1, &F, &X1new, 0, 0, 1, 1);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "FMI");
  dpd_contract222(&F, &X1, &X1new, 1, 1, -1, 1);
  dpd_file2_close(&F);

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "(2WmBeJ + WmBEj) (jb,me)");
  dpd_contract422(&W, &X1, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&W);

  dpd_file2_close(&X1);


  /*** S-D ***/

  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "X_%s_%1s_(2IjAb-IjbA) (%5.3f)", pert, cart, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&F, &X2, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&X2);
  dpd_file2_close(&F);

  sprintf(lbl, "X_%s_%1s_IjAb (%5.3f)", pert, cart, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_contract442(&X2, &W, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
  dpd_contract442(&W, &X2, &X1new, 3, 3, 1, 1);
  dpd_buf4_close(&W);

  dpd_buf4_close(&X2);

  if(params.local && local.filter_singles) local_filter_T1(&X1new, omega);
  else denom1(&X1new, omega);
  dpd_file2_close(&X1new);
}
