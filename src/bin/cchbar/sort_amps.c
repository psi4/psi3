#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void sort_amps(void)
{
  dpdbuf4 t2;

  /* T(iJ,aB) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_sort(&t2, CC_TAMPS, qpsr, 0, 5, "tiJaB");
  dpd_buf4_close(&t2);

  /* TIJAB (IA,JB) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tIAJB");
  dpd_buf4_close(&t2);

  /* Tijab (ia,jb) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tiajb");
  dpd_buf4_close(&t2);

  /* TIjAb (IA,jb) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tIAjb");
  dpd_buf4_close(&t2);

  /* TIjAb (ij,JB) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tiaJB");
  dpd_buf4_close(&t2);

  /* TIjAb (Ib,jA) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_sort(&t2, CC_TAMPS, psrq, 10, 10, "tIbjA");
  dpd_buf4_close(&t2);
  /* TIjAb (jA,Ib) */
  dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tjAIb");
  dpd_buf4_close(&t2);
}
