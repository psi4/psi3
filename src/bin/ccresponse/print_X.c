/*! \file print_X.c
    \ingroup (CCRESPONSE)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void print_X(char *pert, char *cart, int irrep, double omega)
{
  dpdfile2 X1;
  dpdbuf4 X2;
  char lbl[32];

  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert, cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_print(&X1, outfile);
  dpd_file2_close(&X1);
  sprintf(lbl, "X_%s_%1s_IjAb (%5.3f)", pert, cart, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_print(&X2, outfile, 1);
  dpd_buf4_close(&X2);
}
