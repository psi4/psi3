#include <stdio.h>
#include <stdlib.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void e_sort(void)
{
  dpdbuf4 E;

  if(params.ref == 2) {  /** UHF **/
    /*** AA ***/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 21, 0, 21, 0, 0, "E <AI|JK>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 0, 20, "E <IJ|KA>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 2, 20, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_close(&E);

    /*** BB ***/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 31, 10, 31, 10, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 10, 30, "E <ij|ka>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 12, 30, "E <ij||ka> (i>j,ka)");
    dpd_buf4_close(&E);

  }
  else {  /** RHF/ROHF **/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 0, 10, "E <ij|ka>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 2, 10, "E <ij||ka> (i>j,ka)");
    dpd_buf4_close(&E);

    /* <ij|ka> (ij,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 0, 11, "E <ij|ka> (ij,ak)");
    dpd_buf4_close(&E);

    /* <ij||ka> (ij,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 2, 11, "E <ij||ka> (i>j,ak)");
    dpd_buf4_close(&E);

    /* <ia|jk> */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, qpsr, 10, 0, "E <ia|jk>");
    dpd_buf4_close(&E);
  }
}
