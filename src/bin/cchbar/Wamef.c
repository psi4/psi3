#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wamef_build(void) {
  dpdbuf4 Wamef, WAMEF, WAmEf, WaMeF;
  dpdbuf4 F, D_a, D;
  dpdfile2 tia, tIA;
  
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  /* F(ma,e>f) --> W(ma,e>f) */
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_copy(&F, CC_HBAR, "WAMEF");
  dpd_buf4_copy(&F, CC_HBAR, "Wamef");
  dpd_buf4_close(&F);

  /* D(NM,E>F) * T(N,A) --> W(MA,E>F) */
  dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
  dpd_buf4_init(&D_a, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract424(&D_a, &tIA, &WAMEF, 1, 0, 1, 1.0, -1.0);
  dpd_buf4_close(&D_a);
  dpd_buf4_close(&WAMEF);  

  /* D(nm,e>f) * T(n,a) --> W(ma,e>f) */
  dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
  dpd_buf4_init(&D_a, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract424(&D_a, &tia, &Wamef, 1, 0, 1, 1.0, -1.0);
  dpd_buf4_close(&D_a);
  dpd_buf4_close(&Wamef); 

  /* F(ma,ef) --> W(ma,ef) */
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&F, CC_HBAR, pqsr, 10, 5, "WAmEf");
  dpd_buf4_sort(&F, CC_HBAR, pqsr, 10, 5, "WaMeF");
  dpd_buf4_close(&F);

  /* D(ij,ab) --> D(ji,ab) */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_TMP0, qprs, 0, 5, "D <ij|ab> (ji,ab)");
  dpd_buf4_close(&D);

  /* D(mN,Ef) * T(N,A) --> W(mA,Ef) */
  dpd_buf4_init(&WAmEf, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  /* D(Mn,eF) * T(n,a) --> W(Ma,eF) */
  dpd_buf4_init(&WaMeF, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
  dpd_buf4_init(&D, CC_TMP0, 0, 0, 5, 0, 5, 0, "D <ij|ab> (ji,ab)");
  dpd_contract424(&D, &tIA, &WAmEf, 1, 0, 1, -1.0, 1.0);
  dpd_contract424(&D, &tia, &WaMeF, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&WAmEf);
  dpd_buf4_close(&WaMeF);

  dpd_file2_close(&tIA);
  dpd_file2_close(&tia);

  return;
}
