#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmnie_build(void) {
  dpdbuf4 Wmnie, WMNIE, WMnIe, WmNiE, WMniE, WmNIe;
  dpdbuf4 E;
  dpdbuf4 D, D_a;
  dpdfile2 t1;

  /* E(M>N,EI) --> W(M>N,EI) */
  dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_sort(&E, CC_HBAR, pqsr, 2, 11, "WMNIE");
  dpd_buf4_sort(&E, CC_HBAR, pqsr, 2, 11, "Wmnie");
  dpd_buf4_close(&E);


  /* D(M>N,EF) * T(I,F) --> W(M>N,EI) */
  dpd_buf4_init(&WMNIE, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE");
  dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&D_a,&t1,&WMNIE, 3, 1, 0, -1, 1);
  dpd_file2_close(&t1);
  dpd_buf4_close(&D_a);
  dpd_buf4_close(&WMNIE);


  /* D(m>n,ef) * T(i,f) --> W(m>n,ei) */
  dpd_buf4_init(&Wmnie, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie");
  dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&D_a, &t1, &Wmnie, 3, 1, 0, -1, 1);
  dpd_file2_close(&t1);
  dpd_buf4_close(&D_a);
  dpd_buf4_close(&Wmnie);


  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_copy(&E, CC_TMP0, "WMnIe (Mn,Ie)");
  dpd_buf4_copy(&E, CC_TMP1, "WmNiE (mN,iE)");
  dpd_buf4_close(&E);

  /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
  dpd_buf4_init(&WMnIe, CC_TMP0, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&t1, &D, &WMnIe, 1, 2, 1, 1, 1);
  dpd_file2_close(&t1);
  dpd_buf4_close(&D);
  /* W(Mn,Ie) --> W(Mn,eI) */
  dpd_buf4_sort(&WMnIe, CC_HBAR, pqsr, 0, 11, "WMnIe");
  dpd_buf4_close(&WMnIe);

  /* D(mN,fE) * T(i,f) --> W(mN.iE) */
  dpd_buf4_init(&WmNiE, CC_TMP1, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&t1,&D,&WmNiE, 1, 2, 1, 1, 1);
  dpd_file2_close(&t1);
  dpd_buf4_close(&D);
  /* W(mN,iE) --> W(mN,Ei) */
  dpd_buf4_sort(&WmNiE, CC_HBAR, pqsr, 0, 11, "WmNiE");
  dpd_buf4_close(&WmNiE);

  return;
}
