#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void c_sort(void)
{
  dpdbuf4 C, D;

  /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_copy(&C, CC_CINTS, "C <ia||jb>");
  dpd_buf4_close(&C);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_TMP0, psqr, 10, 10, "D <ij|ab> (ib,ja)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_TMP0, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_axpy(&D, &C, -1);
  dpd_buf4_close(&D);
  dpd_buf4_close(&C);

  /* <ia|jb> (bi,ja) */
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort(&C, CC_CINTS, sprq, 11, 10, "C <ia|jb> (bi,ja)");
  dpd_buf4_close(&C);

  /* <ia||jb> (bi,ja) */
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_sort(&C, CC_CINTS, sprq, 11, 10, "C <ia||jb> (bi,ja)");
  dpd_buf4_close(&C);

  /* <ia|jb> (ia,bj) */
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort(&C, CC_CINTS, pqsr, 10, 11, "C <ia|jb> (ia,bj)");
  dpd_buf4_close(&C);

  /* <ia||jb> (ia,bj) */
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_sort(&C, CC_CINTS, pqsr, 10, 11, "C <ia||jb> (ia,bj)");
  dpd_buf4_close(&C);

}
