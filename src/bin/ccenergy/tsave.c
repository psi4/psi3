#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void tsave(void)
{
  struct oe_dpdfile t1;
  struct dpdbuf t2;

  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_copy(&t1, CC_OEI, "tIA", 0, outfile);
  dpd_oe_file_close(&t1);

  dpd_oe_file_init(&t1, CC_OEI, 0, 1, "New tia", 0, outfile);
  dpd_oe_copy(&t1, CC_OEI, "tia", 0, outfile);
  dpd_oe_file_close(&t1);

  dpd_buf_init(&t2, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_copy(&t2, CC_TAMPS, "tIJAB", 0, outfile);
  dpd_buf_close(&t2);

  dpd_buf_init(&t2, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_copy(&t2, CC_TAMPS, "tijab", 0, outfile);
  dpd_buf_close(&t2);

  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_copy(&t2, CC_TAMPS, "tIjAb", 0, outfile);
  dpd_buf_close(&t2);
}
