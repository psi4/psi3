#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void DT2(void)
{
  struct dpdbuf D;

  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0,
	       "D <ij||ab> (i>j,a>b)", 0, outfile);
  dpd_copy(&D, CC_TAMPS, "New tIJAB", 0, outfile);
  dpd_copy(&D, CC_TAMPS, "New tijab", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_TAMPS, "New tIjAb", 0, outfile);
  dpd_buf_close(&D);
}
