#include <stdio.h>
#include <stdlib.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

void BT2(void)
{
  int h;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 B_anti, B;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;
  dpdbuf4 Z1,Z2;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);

    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)");

    dpd_buf4_init(&Z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z2, &newtIjAb, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Z1);

    dpd_buf4_close(&B);

    dpd_buf4_close(&tauIjAb);

    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    dpd_buf4_init(&B_anti, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");

    /* AA and BB terms */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(ab,ij)");

    dpd_contract444(&B_anti, &tauIJAB, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    dpd_buf4_axpy(&Z2, &newtIJAB, 1);
    dpd_buf4_close(&Z2);

    dpd_contract444(&B_anti, &tauijab, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    dpd_buf4_axpy(&Z2, &newtijab, 1);
    dpd_buf4_close(&Z2);

    dpd_buf4_close(&Z1);

    /* AB term */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z2, &newtIjAb, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Z1);

    dpd_buf4_close(&B_anti);  
    dpd_buf4_close(&B);

    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&tauIjAb);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    dpd_contract444(&B, &tauIJAB, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 2, 7, "Z(IJ,AB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(IJ,AB)");
    dpd_buf4_axpy(&Z2, &newtIJAB, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&B);

    dpd_buf4_init(&B, CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 17, 12, 17, 12, 0, "Z(ab,ij)");
    dpd_contract444(&B, &tauijab, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 12, 17, "Z(ij,ab)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 12, 17, 12, 17, 0, "Z(ij,ab)");
    dpd_buf4_axpy(&Z2, &newtijab, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&B);

    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    dpd_contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 22, 28, "Z(Ij,Ab)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z2, &newtIjAb, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&B);

    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&tauIjAb);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }

}
