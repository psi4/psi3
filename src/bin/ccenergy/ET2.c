#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void ET2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 E, t2, t2a, t2b;

/*  timer_on("ET2", outfile); */

  dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
  dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  /*** AA ***/
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
  dpd_buf4_sort(&t2, CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
  dpd_buf4_close(&t2);
  dpd_buf4_init(&t2a, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_buf4_init(&t2b, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
  dpd_buf4_axpy(&t2b, &t2a, -1);
  dpd_buf4_close(&t2b);
  dpd_buf4_axpy(&t2a, &newtIJAB, 1);
  dpd_buf4_close(&t2a);
  dpd_buf4_close(&E);

  /*** BB ***/
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
  dpd_buf4_sort(&t2, CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
  dpd_buf4_close(&t2);
  dpd_buf4_init(&t2a, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_buf4_init(&t2b, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
  dpd_buf4_axpy(&t2b, &t2a, -1);
  dpd_buf4_close(&t2b);
  dpd_buf4_axpy(&t2a, &newtijab, 1);
  dpd_buf4_close(&t2a);
  dpd_buf4_close(&E);

  /*** AB ***/

  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_contract424(&E, &tia, &newtIjAb, 1, 0, 0, -1, 1);
  dpd_buf4_close(&E);
  dpd_buf4_init(&E, CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
  dpd_contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
  dpd_buf4_close(&E);

  dpd_file2_close(&tIA); dpd_file2_close(&tia);

  dpd_buf4_close(&newtIJAB);
  dpd_buf4_close(&newtijab);
  dpd_buf4_close(&newtIjAb);

/*  timer_off("ET2", outfile); */
}
