#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WabeiL1(void)
{
  struct oe_dpdfile newL1;
  struct dpdbuf W, L2;

  dpd_oe_file_init(&newL1, CC_OEI, 0, 1, "New L(I,A)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "W(AM,EF)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "L2(IM,EF)", 0, outfile);
  dpd_contract122(&L2, &W, &newL1, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&W, CC_AbEi, 11, 5, 11, 5, 0, "W(Am,Ef)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "L2(Im,Ef)", 0, outfile);
  dpd_contract122(&L2, &W, &newL1, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_oe_file_close(&newL1);

  dpd_oe_file_init(&newL1, CC_OEI, 0, 1, "New L(i,a)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "W(am,ef)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "L2(im,ef)", 0, outfile);
  dpd_contract122(&L2, &W, &newL1, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&W, CC_aBeI, 11, 5, 11, 5, 0, "W(aM,eF)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "L2(iM,eF)", 0, outfile);
  dpd_contract122(&L2, &W, &newL1, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_oe_file_close(&newL1);
}
