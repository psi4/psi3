#include <stdio.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void d_sort(void)
{
  struct dpdbuf D;

  dpd_buf_init(&D, CC_DINTS, 2, 7, 0, 5, 1, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_DINTS, "D <ij||ab> (i>j,a>b)", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 2, 5, 0, 5, 1, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_DINTS, "D <ij||ab> (i>j,ab)", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 7, 0, 5, 1, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_DINTS, "D <ij||ab> (ij,a>b)", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 1, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_DINTS, "D <ij||ab>", 0, outfile);
  dpd_buf_close(&D);


  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap23(&D, CC_DINTS, 10, 10, "D <ij|ab> (ia,jb)", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_swap12(&D, CC_DINTS, 11, 10, "D <ij|ab> (ai,jb)", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij||ab>", 0, outfile);
  dpd_swap23(&D, CC_DINTS, 10, 10, "D <ij||ab> (ia,jb)", 0, outfile);
  dpd_buf_close(&D);
  
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_swap24(&D, CC_DINTS, 10, 10, "D <ij|ab> (ib,ja)",  0, outfile);
  dpd_buf_close(&D);
}
