#include <stdio.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void d_sort(void)
{
  dpdbuf4 D;

  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 0, 5, 1, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,ab)");
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 5, 1, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (ij,a>b)");
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 1, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab>");
  dpd_buf4_close(&D);

  /* <ij|ab> (ia,jb) */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_DINTS, prqs, 10, 10, "D <ij|ab> (ia,jb)");
  dpd_buf4_close(&D);

  /* <ij|ab> (ai,jb) */
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_sort(&D, CC_DINTS, qprs, 11, 10, "D <ij|ab> (ai,jb)");
  dpd_buf4_close(&D);

  /* <ij||ab> (ia,jb) */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_buf4_sort(&D, CC_DINTS, prqs, 10, 10, "D <ij||ab> (ia,jb)");
  dpd_buf4_close(&D);
  
  /* <ij|ab> (ib,ja) */
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_sort(&D, CC_DINTS, psrq, 10, 10, "D <ij|ab> (ib,ja)");
  dpd_buf4_close(&D);

  /* <ij|ab> (ib,aj) */
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ib,aj)");
  dpd_buf4_close(&D);

  /* <ij|ab> (ia,bj) */
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ia,bj)");
  dpd_buf4_close(&D);

  /* <ij||ab> (ia,bj) */
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij||ab> (ia,bj)");
  dpd_buf4_close(&D);
}
