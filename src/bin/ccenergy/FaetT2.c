#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FaetT2(void)
{
  dpdfile2 FAEt, Faet;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2;

  dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
  dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

  dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
  dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");

  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&tIJAB, &FAEt, &t2, 3, 1, 0, 1, 0);
  dpd_contract244(&FAEt, &tIJAB, &t2, 1, 2, 1, 1, 1);
  dpd_buf4_axpy(&t2, &newtIJAB, 1);
  dpd_buf4_close(&t2);

  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&tijab, &Faet, &t2, 3, 1, 0, 1, 0);
  dpd_contract244(&Faet, &tijab, &t2, 1, 2, 1, 1, 1);
  dpd_buf4_axpy(&t2, &newtijab, 1);
  dpd_buf4_close(&t2);

  dpd_contract424(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1);
  dpd_contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);

  dpd_file2_close(&FAEt);  dpd_file2_close(&Faet);

  dpd_buf4_close(&tIJAB);
  dpd_buf4_close(&tijab);
  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&newtIJAB);
  dpd_buf4_close(&newtijab);
  dpd_buf4_close(&newtIjAb);
}
