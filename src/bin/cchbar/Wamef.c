#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* Wamef_build(): Computes all contributions to the Wamef HBAR matrix
** elements, whose spin-orbital definition is:
**
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** The storage and naming convention for each spin case are 
** as follows:
**
** Spin Case     Storage      Name
** ----------    ---------    --------
** WAMEF         (MA,E>F)      "WAMEF"
** Wamef         (ma,e>f)      "Wamef"
** WAmEf         (mA,Ef)       "WAmEf"
** WaMeF         (Ma,eF)       "WaMeF"
** -----------------------------------
**
** TDC, June 2002
**
** RHF Cases:  Note that only the WAmEf spin case is required, and
** we store it AS WRITTEN, (Am,Ef).
**
** TDC, March 2004
*/

void Wamef_build(void) {
  dpdbuf4 Wamef, WAMEF, WAmEf, WaMeF, W;
  dpdbuf4 F, D_a, D;
  dpdfile2 tia, tIA;

  if(params.ref == 0) { 

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_copy(&F, CC_HBAR, "WAmEf");
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&W);

  }
  else if(params.ref == 1) { /** ROHF **/
  
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

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /* <AM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_buf4_copy(&F, CC_HBAR, "WAMEF");
    dpd_buf4_close(&F);

    /* <am||ef> --> W(am,e>f) */
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC_HBAR, "Wamef");
    dpd_buf4_close(&F);

    /* T(N,A) <NM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&WAMEF, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract244(&tIA, &D, &WAMEF, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAMEF);  

    /* T(n,a) <nm||ef> --> W(am,e>f) */
    dpd_buf4_init(&Wamef, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tia, &D, &Wamef, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Wamef); 

    /* <Am|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_buf4_copy(&F, CC_HBAR, "WAmEf");
    dpd_buf4_close(&F);

    /* <aM|eF> --> W(aM,eF) */
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_buf4_copy(&F, CC_HBAR, "WaMeF");
    dpd_buf4_close(&F);

    /* T(N,A) <Nm|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&WAmEf, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract244(&tIA, &D, &WAmEf, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAmEf);

    /* T(n,a) <nM|eF> --> W(aM,eF) */
    dpd_buf4_init(&WaMeF, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract244(&tia, &D, &WaMeF, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WaMeF);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  } /** UHF **/

  return;
}
