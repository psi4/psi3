#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

/* Wmbej_build(): Build the Wmbej intermediate.
** 
** Wmbej = <mb||ej> + t(j,f) <mb||ef> - t(n,b) <mn||ej>
**         - [1/2 t(jn,fb) + t(j,f) t(n,b)] <mn||ef>
**
** Spin cases for UHF and ROHF orbitals:
** -------------------------------------
**
**
** TDC
** May 2000
*/

void Wmbej_build(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ; 
  dpdbuf4 C, D, E, F, X, tIAjb, tiaJB, t2, W, Y;
  dpdfile2 tIA, tia;

  timer_on("C->Wmbej");

  /* W(mb,je) <-- <mb||ej> */

  if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <Ia|Jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    }
  else { /*** RHF/ROHF ***/


    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMBEJ", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "WmBEj", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "WMbEj");
    dpd_buf4_copy(&D, CC_TMP0, "WmBeJ");
    dpd_buf4_close(&D);

  }

  timer_off("C->Wmbej");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  timer_on("F->Wmbej");
  
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 7, 0, "F <ia||bc> (ia,b>c)");
  dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
  dpd_contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
  dpd_buf4_close(&WMBEJ);
  dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
  dpd_contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
  dpd_buf4_close(&Wmbej);
  dpd_buf4_close(&F);

  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
  dpd_contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
  dpd_buf4_close(&WMbEj);
  dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
  dpd_contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
  dpd_buf4_close(&WmBeJ);

  dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
  dpd_buf4_close(&WMbeJ);
  dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
  dpd_buf4_close(&WmBEj);
  dpd_buf4_close(&F);

  timer_off("F->Wmbej");

  timer_on("E->Wmbej");

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
  dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
  dpd_contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
  dpd_buf4_close(&WMBEJ);
  dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
  dpd_contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
  dpd_buf4_close(&Wmbej);
  dpd_buf4_close(&E);

  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
  dpd_contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
  dpd_buf4_close(&WMbEj);
  dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
  dpd_contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
  dpd_buf4_close(&WmBeJ);
  dpd_buf4_close(&E);

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
  dpd_buf4_close(&WMbeJ);
  dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
  dpd_buf4_close(&WmBEj);
  dpd_buf4_close(&E);

  timer_off("E->Wmbej");

  dpd_file2_close(&tIA);  dpd_file2_close(&tia);


  /* Convert to (ME,JB) for remaining terms */

  timer_on("sort Wmbej");

  dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
  dpd_buf4_sort(&WMBEJ, CC_HBAR, prsq, 10, 10, "WMBEJ");
  dpd_buf4_close(&WMBEJ);

  dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
  dpd_buf4_sort(&Wmbej, CC_HBAR, prsq, 10, 10, "Wmbej");
  dpd_buf4_close(&Wmbej);

  dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
  dpd_buf4_sort(&WMbEj, CC_HBAR, prsq, 10, 10, "WMbEj");
  dpd_buf4_close(&WMbEj);

  dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
  dpd_buf4_sort(&WmBeJ, CC_HBAR, prsq, 10, 10, "WmBeJ");
  dpd_buf4_close(&WmBeJ);

  dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_buf4_sort(&WMbeJ, CC_HBAR, psrq, 10, 10, "WMbeJ");
  dpd_buf4_close(&WMbeJ);

  dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_buf4_sort(&WmBEj, CC_HBAR, psrq, 10, 10, "WmBEj");
  dpd_buf4_close(&WmBEj);

  timer_off("sort Wmbej");


  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  /*** AAAA ***/

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
  dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  /*** BBBB ***/

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
  dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
  dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  /*** ABAB ***/
  
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
  dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  /*** BABA ***/

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&t2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
  dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  /*** ABBA ***/
  
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&D);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
  dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  /*** BAAB ***/

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  timer_on("X->Wmbej");
  dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
  timer_off("X->Wmbej");
  dpd_buf4_close(&D);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
  dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y);

  dpd_file2_close(&tIA);
  dpd_file2_close(&tia);

}
