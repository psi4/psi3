#include <stdio.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void e_sort(void)
{
  dpdbuf4 E;

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
