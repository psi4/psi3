#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wamef_build(void) {
  struct dpdbuf Wamef, WAMEF, WAmEf, WaMeF;
  struct dpdbuf F, D_a, D;
  struct oe_dpdfile tia, tIA;
  
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  /* F(ma,e>f) --> W(ma,e>f) */
  dpd_buf_init(&F, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_copy(&F, CC_HBAR, "WAMEF", 0, outfile);
  dpd_copy(&F, CC_HBAR, "Wamef", 0, outfile);
  dpd_buf_close(&F);

  /* D(NM,E>F) * T(N,A) --> W(MA,E>F) */
  dpd_buf_init(&WAMEF, CC_HBAR, 10, 7, 10, 7, 0, "WAMEF", 
               0, outfile);
  dpd_buf_init(&D_a, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)", 
               0, outfile);
  dpd_contract221(&D_a, &tIA, &WAMEF, 1, 0, 1, 1.0, -1.0, 0, outfile);
  dpd_buf_close(&D_a);
  dpd_buf_close(&WAMEF);  

  /* D(nm,e>f) * T(n,a) --> W(ma,e>f) */
  dpd_buf_init(&Wamef, CC_HBAR, 10, 7, 10, 7, 0, "Wamef", 0, outfile);
  dpd_buf_init(&D_a, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)", 
               0, outfile);
  dpd_contract221(&D_a, &tia, &Wamef, 1, 0, 1, 1.0, -1.0, 0, outfile);
  dpd_buf_close(&D_a);
  dpd_buf_close(&Wamef); 

  /* F(ma,ef) --> W(ma,ef) */
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap34(&F, CC_HBAR, 10, 5, "WAmEf", 0, outfile);
  dpd_swap34(&F, CC_HBAR, 10, 5, "WaMeF", 0, outfile);
  dpd_buf_close(&F);

  /* D(ij,ab) --> D(ji,ab) */
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap12(&D, CC_TMP0, 0, 5, "D <ij|ab> (ji,ab)", 0, outfile);
  dpd_buf_close(&D);

  /* D(mN,Ef) * T(N,A) --> W(mA,Ef) */
  dpd_buf_init(&WAmEf, CC_HBAR, 10, 5, 10, 5, 0, "WAmEf", 0, outfile);
  /* D(Mn,eF) * T(n,a) --> W(Ma,eF) */
  dpd_buf_init(&WaMeF, CC_HBAR, 10, 5, 10, 5, 0, "WaMeF", 0, outfile);
  dpd_buf_init(&D, CC_TMP0, 0, 5, 0, 5, 0, "D <ij|ab> (ji,ab)", 0, outfile);
  dpd_contract221(&D, &tIA, &WAmEf, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_contract221(&D, &tia, &WaMeF, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&WAmEf);
  dpd_buf_close(&WaMeF);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);

  return;
}
