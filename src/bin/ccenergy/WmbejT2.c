#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void WmbejT2(void)
{
  struct dpdbuf T2new, T2, W;

/*  timer_on("WmbejT2", outfile); */


  /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (IA,JB) Temp",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  /* P(IJ) T2(IA,JB) */
  dpd_swap13(&T2new, CC_TMP1, 10, 10, "T2 (JA,IB) Temp", 0, outfile);
  /* P(AB) T2(IA,JB) */
  dpd_swap24(&T2new, CC_TMP2, 10, 10, "T2 (IB,JA) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  /* P(IJ) P(AB) T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (JA,IB) Temp",
	       0, outfile);
  dpd_swap24(&T2new, CC_TMP3, 10, 10, "T2 (JB,IA) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (IA,JB) Temp",
	       0, outfile);
  dpd_buf_init(&T2, CC_TMP1, 10, 10, 10, 10, 0, "T2 (JA,IB) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP2, 10, 10, 10, 10, 0, "T2 (IB,JA) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP3, 10, 10, 10, 10, 0, "T2 (JB,IA) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, 1, 0, outfile);
  dpd_buf_close(&T2);
  /* T2(IA,JB) --> T2(IJ,AB) */
  dpd_swap23(&T2new, CC_TMP1, 0, 5, "T2 (IJ,AB) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (IJ,AB) Temp", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);


  
  /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (ia,jb) Temp",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  /* P(ij) T2(ia,jb) */
  dpd_swap13(&T2new, CC_TMP1, 10, 10, "T2 (ja,ib) Temp", 0, outfile);
  /* P(ab) T2(ia,jb) */
  dpd_swap24(&T2new, CC_TMP2, 10, 10, "T2 (ib,ja) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  /* P(ij) P(ab) T2(ia,jb) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (ja,ib) Temp",
	       0, outfile);
  dpd_swap24(&T2new, CC_TMP3, 10, 10, "T2 (jb,ia) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  /* T2(ia,jb) - T2(ja,ib) - T2(ib,ja) + T2(jb,ia) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (ia,jb) Temp",
	       0, outfile);
  dpd_buf_init(&T2, CC_TMP1, 10, 10, 10, 10, 0, "T2 (ja,ib) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP2, 10, 10, 10, 10, 0, "T2 (ib,ja) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP3, 10, 10, 10, 10, 0, "T2 (jb,ia) Temp", 0, outfile);
  dpd_axpy(&T2, &T2new, 1, 0, outfile);
  dpd_buf_close(&T2);
  /* T2(ia,jb) --> T2(ij,ab) */
  dpd_swap23(&T2new, CC_TMP1, 0, 5, "T2 (ij,ab) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (ij,ab) Temp", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "New tijab", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);


  

  /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (IA,jb)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract222(&W, &T2, &T2new, 1, 0, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_contract222(&W, &T2, &T2new, 1, 0, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* T2(IA,jb) --> T2(Ij,Ab) */
  dpd_swap23(&T2new, CC_TMP1, 0, 5, "T2 (Ij,Ab) Temp", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (Ij,Ab) Temp", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);

  /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
  dpd_buf_init(&T2new, CC_TMP0, 10, 10, 10, 10, 0, "T2 (Ib,jA)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_contract222(&T2, &W, &T2new, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_contract222(&W, &T2, &T2new, 1, 0, 1, 1, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);

  /* T2(Ib,jA) --> T2(Ij,bA) --> T2(Ij,Ab) */
  dpd_swap23(&T2new, CC_TMP1, 0, 5, "T2 (Ij,bA)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (Ij,bA)", 0, outfile);
  dpd_swap34(&T2new, CC_TMP0, 0, 5, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);

/*  timer_off("WmbejT2", outfile); */
}
