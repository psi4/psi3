#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmbej_build(void) {
  struct dpdbuf WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  struct dpdbuf tIAJB, tjAIb, tiajb, tIAjb, tiaJB, tIbjA;
  struct dpdbuf D;

  /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 
               0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)", 
               0, outfile);
  dpd_buf_init(&tiajb, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_contract222(&D, &tiajb, &Wmbej, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiajb);

  /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
  dpd_buf_init(&tiaJB, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiaJB);

  dpd_buf_close(&Wmbej);


  /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 
               0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)", 
               0, outfile);
  dpd_buf_init(&tIAJB, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_contract222(&D, &tIAJB, &WMBEJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAJB);

  /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
  dpd_buf_init(&tIAjb, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAjb);

  dpd_buf_close(&WMBEJ);


  /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 
               0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)", 
               0, outfile);
  dpd_buf_init(&tIAjb, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAjb);

  /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&tIAJB, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_contract222(&D, &tIAJB, &WmBeJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAJB);

  dpd_buf_close(&WmBeJ);


  /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
  dpd_buf_init(&WMbEj, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)", 
               0, outfile);
  dpd_buf_init(&tiaJB, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiaJB);

  /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&tiajb, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_contract222(&D, &tiajb, &WMbEj, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiajb);

  dpd_buf_close(&WMbEj);


  /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)",
	       0, outfile);
  dpd_buf_init(&tjAIb, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_contract222(&D, &tjAIb, &WmBEj, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tjAIb);

  dpd_buf_close(&WmBEj);


  /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)",
	       0, outfile);
  dpd_buf_init(&tIbjA, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_contract222(&D, &tIbjA, &WMbeJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIbjA);

  dpd_buf_close(&WMbeJ);

  return;
}

