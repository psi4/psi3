#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

/* WmbejT2(): Contributions of Wmbej intermediates to T2.
**
** t(ij,ab) <--- P(ij) P(ab) t(im,ae) Wmbej
**
** Spin cases for UHF or ROHF orbitals:
** ------------------------------------
**                 *** AA ***
** + T(IA,ME) W(ME,JB) + T(IA,me) W(me,JB)
** - T(JA,ME) W(ME,IB) - T(JA,me) W(me,IB)
** - T(IB,ME) W(ME,JA) - T(IB,me) W(me,JA)
** + T(JB,ME) W(ME,IA) + T(JB,me) W(me,IA)
**
**                 *** BB ***
** + T(ia,me) W(me,jb) + T(ia,ME) W(ME,jb)
** - T(ja,me) W(me,ib) - T(ja,ME) W(ME,ib)
** - T(ib,me) W(me,ja) - T(ib,ME) W(ME,ja)
** + T(jb,me) W(me,ia) + T(jb,ME) W(ME,ia)
**
**                 *** AB ***
** + T(IA,ME) W(ME,jb) + T(IA,me) W(me,jb)
** + T(MA,je) W(Me,Ib) + T(IE,mb) W(mE,jA)
** + T(jb,ME) W(ME,IA) + T(jb,me) W(me,IA)
**
** For the AA and BB spin cases, only the first two terms of each need to
** be evaluated, while for AB, all six terms are different.
** 
** The current version of this code requires ten contractions (two each for
** AA and BB and six for AB), one complex sort each for AA and BB, two
** complex sorts for AB, and three simple sorts each for AA and BB.
**
** For RHF orbitals:
** -----------------
**
**
** TDC
** May 2000
** Revised August 2001
*/

void WmbejT2(void)
{
  dpdbuf4 T2new, T2, W, T2B, W1, W2, Z;

  if(params.ref == 0) { /** RHF **/
    /*** AB ***/

    /* 2 W(ME,jb) + W(Me,Jb) */
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_copy(&W, CC_TMP0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W1, CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_axpy(&W2, &W1, 2);
    dpd_buf4_close(&W2);
    dpd_buf4_close(&W1);


    /* T2(Ib,mE) * W(mE,jA) --> Z(Ib,jA) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    /* T2(Ib,jA) --> T2(IA,jb) (part III) */
    dpd_buf4_sort(&T2new, CC_TMP0, psrq, 10, 10, "T2 (IA,jb) 3");
    dpd_buf4_close(&T2new);

    /* 1/2 [ (2 T2(IA,me) - T2(IE,ma)) * (2 W(ME,jb) + W(Me,Jb)] --> T2(IA,jb) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 0.5, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    /* 1/2 Z(Ib,jA) + T2(IA,jb) --> T2(IA,jb) (Part I) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    dpd_buf4_axpy(&Z, &T2new, 0.5);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2new);

    /* T2(IA,jb) (I) + T2(IA,jb) (III) --> T2(IA,jb) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    dpd_buf4_init(&T2, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 3");
    dpd_buf4_axpy(&T2, &T2new, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_sort(&T2new, CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) (1+3)");
    dpd_buf4_close(&T2new);

    dpd_buf4_init(&T2new, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    dpd_buf4_sort(&T2new, CC_TMP0, qpsr, 0, 5, "T2 (Ij,Ab) (2+4)");
    dpd_buf4_close(&T2new);

    /* T2(Ij,Ab) <--- I + II + III + IV */
    dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    dpd_buf4_axpy(&T2, &T2new, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (2+4)");
    dpd_buf4_axpy(&T2, &T2new, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_close(&T2new);

  }
  else if(params.ref == 1) { /** ROHF **/

    /*** AA ***/

    /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(IA,JB) --> T2(IJ,AB) */
    dpd_buf4_sort(&T2new, CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    dpd_buf4_close(&T2new);
  
    /* P(IJ) P(AB) T2(IA,JB) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    dpd_buf4_sort(&T2new, CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    dpd_buf4_sort(&T2new, CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    dpd_buf4_sort(&T2new, CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    dpd_buf4_axpy(&T2, &T2new, -1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    dpd_buf4_axpy(&T2, &T2new, -1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    dpd_buf4_axpy(&T2, &T2new, +1);
    dpd_buf4_close(&T2);

    /* T2(IJ,AB) --> New T2(IJ,AB) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&T2new, &T2, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&T2new);


    /*** BB ***/

    /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(ia,jb) --> T2(ij,ab) */
    dpd_buf4_sort(&T2new, CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    dpd_buf4_close(&T2new);
  
    /* P(ij) P(ab) T2(ia,jb) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    dpd_buf4_sort(&T2new, CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    dpd_buf4_sort(&T2new, CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    dpd_buf4_sort(&T2new, CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) + T2(ji,ba) --> T2(ij,ab) */
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    dpd_buf4_axpy(&T2, &T2new, -1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    dpd_buf4_axpy(&T2, &T2new, -1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    dpd_buf4_axpy(&T2, &T2new, +1);
    dpd_buf4_close(&T2);
  
    /* T2(ij,ab) --> New T2(ij,ab) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&T2new, &T2, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&T2new);


    /*** AB ***/

    /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(IA,jb) --> T2(Ij,Ab) (part 1) */
    dpd_buf4_sort(&T2new, CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) 1");
    dpd_buf4_close(&T2new);

    /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
    dpd_buf4_init(&T2new, CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (Ib,jA)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    timer_on("WmbejT2 444");
    dpd_contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    timer_on("WmbejT2 444");
    dpd_contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);

    /* T2(Ib,jA) --> T2(Ij,Ab) (part 2) */
    dpd_buf4_sort(&T2new, CC_TMP0, prsq, 0, 5, "T2 (Ij,Ab) 2");
    dpd_buf4_close(&T2new);


    /* T2(Ij,Ab) (part 1) + T2(Ij,Ab) (part 2) --> New T2(Ij,Ab) */
    dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 1");
    dpd_buf4_axpy(&T2, &T2new, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 2");
    dpd_buf4_axpy(&T2, &T2new, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_close(&T2new);

  }
}
