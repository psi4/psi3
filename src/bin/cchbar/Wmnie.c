#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmnie_build(void) {
  struct dpdbuf Wmnie, WMNIE, WMnIe, WmNiE, WMniE, WmNIe;
  struct dpdbuf E;
  struct dpdbuf D, D_a;
  struct oe_dpdfile t1;

  /* E(M>N,EI) --> W(M>N,EI) */
  dpd_buf_init(&E, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)", 
               0, outfile);
  dpd_swap34(&E, CC_HBAR, 2, 11, "WMNIE", 0, outfile);
  dpd_swap34(&E, CC_HBAR, 2, 11, "Wmnie", 0, outfile);
  dpd_buf_close(&E);


  /* D(M>N,EF) * T(I,F) --> W(M>N,EI) */
  dpd_buf_init(&WMNIE, CC_HBAR, 2, 11, 2, 11, 0, "WMNIE", 0, outfile);
  dpd_buf_init(&D_a, CC_DINTS, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)", 
               0, outfile);
  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&D_a,&t1,&WMNIE, 3, 1, 0, -1, 1, 0, outfile);
  dpd_oe_file_close(&t1);
  dpd_buf_close(&D_a);
  dpd_buf_close(&WMNIE);


  /* D(m>n,ef) * T(i,f) --> W(m>n,ei) */
  dpd_buf_init(&Wmnie, CC_HBAR, 2, 11, 2, 11, 0, "Wmnie", 0, outfile);
  dpd_buf_init(&D_a, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)", 
               0, outfile);
  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&D_a, &t1, &Wmnie, 3, 1, 0, -1, 1, 0, outfile);
  dpd_oe_file_close(&t1);
  dpd_buf_close(&D_a);
  dpd_buf_close(&Wmnie);


  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_copy(&E, CC_TMP0, "WMnIe (Mn,Ie)", 0, outfile);
  dpd_copy(&E, CC_TMP1, "WmNiE (mN,iE)", 0, outfile);
  dpd_buf_close(&E);

  /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
  dpd_buf_init(&WMnIe, CC_TMP0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&t1, &D, &WMnIe, 1, 2, 1, 1, 1, 0, outfile);
  dpd_oe_file_close(&t1);
  dpd_buf_close(&D);
  /* W(Mn,Ie) --> W(Mn,eI) */
  dpd_swap34(&WMnIe, CC_HBAR, 0, 11, "WMnIe", 0, outfile);
  dpd_buf_close(&WMnIe);

  /* D(mN,fE) * T(i,f) --> W(mN.iE) */
  dpd_buf_init(&WmNiE, CC_TMP1, 0, 10, 0, 10, 0, "WmNiE (mN,iE)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&t1,&D,&WmNiE, 1, 2, 1, 1, 1, 0, outfile);
  dpd_oe_file_close(&t1);
  dpd_buf_close(&D);
  /* W(mN,iE) --> W(mN,Ei) */
  dpd_swap34(&WmNiE, CC_HBAR, 0, 11, "WmNiE", 0, outfile);
  dpd_buf_close(&WmNiE);

  return;
}
