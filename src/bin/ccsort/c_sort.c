#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void c_sort(void)
{
  struct dpdbuf C, D;

  /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_copy(&C, CC_CINTS, "C <ia||jb>", 0, outfile);
  dpd_buf_close(&C);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_sort(&D, CC_TMP0, pqsr, 0, 5, "D <ij|ab> (ij,ba)", 0, outfile);
/*  dpd_swap34(&D, CC_TMP0, 0, 5, "D <ij|ab> (ij,ba)", 0, outfile); */
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_TMP0, 0, 5, 0, 5, 0, "D <ij|ab> (ij,ba)", 0, outfile);
  dpd_buf_sort(&D, CC_TMP1, prqs, 10, 10, "D <ij|ab> (ib,ja)", 0, outfile);
/*  dpd_swap23(&D, CC_TMP1, 10, 10, "D <ij|ab> (ib,ja)", 0, outfile); */
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_TMP1, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_axpy(&D, &C, -1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&C);

  /* <ia|jb> (bi,ja) */ 
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_swap24(&C, CC_TMP0, 10, 10, "C <ia|jb> (ib,ja)", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&C, CC_TMP0, 10, 10, 10, 10, 0, "C <ia|jb> (ib,ja)", 0, outfile);
  dpd_buf_sort(&C, CC_CINTS, qprs, 11, 10, "C <ia|jb> (bi,ja)", 0, outfile);
/*   dpd_swap12(&C, CC_CINTS, 11, 10, "C <ia|jb> (bi,ja)", 0, outfile); */
  dpd_buf_close(&C);

  /* <ia||jb> (bi,ja) */
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>",
	       0, outfile);
  dpd_swap24(&C, CC_TMP0, 10, 10, "C <ia||jb> (ib,ja)", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&C, CC_TMP0, 10, 10, 10, 10, 0, "C <ia||jb> (ib,ja)",
	       0, outfile);
  dpd_swap12(&C, CC_CINTS, 11, 10, "C <ia||jb> (bi,ja)", 0, outfile);
  dpd_buf_close(&C);
}
