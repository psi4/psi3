#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Wamef_build(void) {
  dpdbuf4 Wamef, WAMEF, WAmEf, WaMeF;
  dpdbuf4 F, D_a, D;
  dpdfile2 tia, tIA;

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
  
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

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /* F(ma,e>f) --> W(ma,e>f) */
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_buf4_copy(&F, CC_HBAR, "WAMEF");
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_buf4_copy(&F, CC_HBAR, "Wamef");
    dpd_buf4_close(&F);

    /* D(NM,E>F) * T(N,A) --> W(MA,E>F) */
    dpd_buf4_init(&WAMEF, CC_HBAR, 0, 20, 7, 20, 7, 0, "WAMEF");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract424(&D_a, &tIA, &WAMEF, 1, 0, 1, 1.0, -1.0);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&WAMEF);  

    /* D(nm,e>f) * T(n,a) --> W(ma,e>f) */
    dpd_buf4_init(&Wamef, CC_HBAR, 0, 30, 17, 30, 17, 0, "Wamef");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract424(&D_a, &tia, &Wamef, 1, 0, 1, 1.0, -1.0);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&Wamef); 

    /* F(ma,ef) --> W(ma,ef) */
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_buf4_sort(&F, CC_HBAR, pqsr, 27, 28, "WAmEf");
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_buf4_sort(&F, CC_HBAR, pqsr, 24, 29, "WaMeF");
    dpd_buf4_close(&F);

    /* D(Ij,Ab) --> D(jI,Ab) */
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_sort(&D, CC_TMP0, qprs, 23, 28, "D <Ij|Ab> (jI,Ab)");
    dpd_buf4_close(&D);

    /* D(iJ,aB) --> D(Ji,aB) */
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_buf4_sort(&D, CC_TMP0, qprs, 22, 29, "D <iJ|aB> (Ji,aB)");
    dpd_buf4_close(&D);

    /* D(mN,Ef) * T(N,A) --> W(mA,Ef) */
    dpd_buf4_init(&WAmEf, CC_HBAR, 0, 27, 28, 27, 28, 0, "WAmEf");
    dpd_buf4_init(&D, CC_TMP0, 0, 23, 28, 23, 28, 0, "D <Ij|Ab> (jI,Ab)");
    dpd_contract424(&D, &tIA, &WAmEf, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAmEf);

    /* D(Mn,eF) * T(n,a) --> W(Ma,eF) */
    dpd_buf4_init(&WaMeF, CC_HBAR, 0, 24, 29, 24, 29, 0, "WaMeF");
    dpd_buf4_init(&D, CC_TMP0, 0, 22, 29, 22, 29, 0, "D <iJ|aB> (Ji,aB)");
    dpd_contract424(&D, &tia, &WaMeF, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WaMeF);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  } /** UHF **/

  return;
}
