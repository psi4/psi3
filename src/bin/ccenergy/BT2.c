#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void BT2(void)
{
  int h;
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf B_anti, B;
  struct dpdbuf tauIJAB, tauijab, tauIjAb;
  struct dpdbuf Z1,Z2;

/*  timer_on("BT2", outfile); */

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);

  dpd_buf_init(&B_anti, CC_BINTS, 7, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);

  /* AA and BB terms */
  dpd_buf_init(&Z1, CC_TMP0, 7, 2, 7, 2, 0, "Z(ab,ij)", 0, outfile);
  dpd_contract222(&B_anti, &tauIJAB, &Z1, 0, 0, 1, 0, 0, outfile);
  dpd_buf_sort(&Z1, CC_TMP0, rspq, 2, 7, "Z(ij,ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 2, 7, 2, 7, 0, "Z(ij,ab)", 0, outfile);
  dpd_axpy(&Z2, &newtIJAB, 1, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_contract222(&B_anti, &tauijab, &Z1, 0, 0, 1, 0, 0, outfile);
  dpd_buf_sort(&Z1, CC_TMP0, rspq, 2, 7, "Z(ij,ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 2, 7, 2, 7, 0, "Z(ij,ab)", 0, outfile);
  dpd_axpy(&Z2, &newtijab, 1, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&Z1);

  /* AB term */
  dpd_buf_init(&Z1, CC_TMP0, 5, 0, 5, 0, 0, "Z(Ab,Ij)", 0, outfile);
  dpd_contract222(&B, &tauIjAb, &Z1, 0, 0, 1, 0, 0, outfile);
  dpd_buf_sort(&Z1, CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 0, 5, 0, 5, 0, "Z(Ij,Ab)", 0, outfile);
  dpd_axpy(&Z2, &newtIjAb, 1, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&Z1);

  dpd_buf_close(&B_anti);  dpd_buf_close(&B);

  dpd_buf_close(&tauIJAB); dpd_buf_close(&tauijab); dpd_buf_close(&tauIjAb);

  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);

/*  timer_off("BT2", outfile); */
}
