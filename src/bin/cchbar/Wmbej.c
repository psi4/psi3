#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmbej_build(void) {
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  dpdbuf4 tIAJB, tjAIb, tiajb, tIAjb, tiaJB, tIbjA;
  dpdbuf4 D;

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    dpd_buf4_init(&Wmbej, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    dpd_buf4_close(&Wmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    dpd_buf4_close(&WMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    dpd_buf4_close(&WmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    dpd_buf4_init(&WMbEj, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    dpd_buf4_close(&WMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    dpd_buf4_init(&WmBEj, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tjAIb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_contract444(&D, &tjAIb, &WmBEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tjAIb);

    dpd_buf4_close(&WmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    dpd_buf4_init(&WMbeJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tIbjA, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&D, &tIbjA, &WMbeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIbjA);

    dpd_buf4_close(&WMbeJ);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    dpd_buf4_init(&Wmbej, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_contract444(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    dpd_buf4_close(&Wmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_contract444(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    dpd_buf4_close(&WMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    dpd_buf4_close(&WmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    dpd_buf4_init(&WMbEj, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    dpd_buf4_close(&WMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    dpd_buf4_init(&WmBEj, CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    dpd_buf4_init(&tjAIb, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    dpd_contract444(&D, &tjAIb, &WmBEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tjAIb);

    dpd_buf4_close(&WmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    dpd_buf4_init(&WMbeJ, CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    dpd_buf4_init(&tIbjA, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_contract444(&D, &tIbjA, &WMbeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIbjA);

    dpd_buf4_close(&WMbeJ);

  } /** UHF **/

  return;
}

