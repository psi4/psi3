/*! \file 
    \ingroup (CCRESPONSE)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#define EXTERN
#include "globals.h"

void denom1(dpdfile2 *X1, double omega);
void denom2(dpdbuf4 *X2, double omega);
void local_filter_T1(dpdfile2 *T1, double omega);
void local_filter_T2(dpdbuf4 *T2, double omega);

void init_X(char *pert, char *cart, int irrep, double omega)
{
  char lbl[32];
  dpdfile2 mu1, X1, FAE, FMI;
  dpdbuf4 X2, mu2;

  sprintf(lbl, "%sBAR_%1s_IA", pert, cart);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert, cart, omega);
  if(!params.restart || !psio_tocscan(CC_OEI, lbl)) {
    dpd_file2_copy(&mu1, CC_OEI, lbl);
    dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
    if(params.local && local.filter_singles) local_filter_T1(&X1, omega);
    else denom1(&X1, omega);
    dpd_file2_close(&X1);
  }
  else fprintf(outfile, "\tUsing existing %s amplitudes.\n", lbl);
  dpd_file2_close(&mu1);

  sprintf(lbl, "%sBAR_%1s_IjAb", pert, cart);
  dpd_buf4_init(&mu2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_%1s_IjAb (%5.3f)", pert, cart, omega);
  if(!params.restart || !psio_tocscan(CC_LR, lbl)) {
    dpd_buf4_copy(&mu2, CC_LR, lbl);
    dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    if(params.local) local_filter_T2(&X2, omega);
    else denom2(&X2, omega);
    dpd_buf4_close(&X2);
  }
  else fprintf(outfile, "\tUsing existing %s amplitudes.\n", lbl);
  dpd_buf4_close(&mu2);
}
