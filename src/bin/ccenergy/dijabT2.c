#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void dijabT2(void)
{
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf dIJAB, dijab, dIjAb;

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&dIJAB, CC_DENOM, 1, 6, 1, 6, 0, "dIJAB", 0, outfile);
  dpd_dirprd(&dIJAB, &newtIJAB, 0, outfile);
  dpd_buf_close(&newtIJAB);
  dpd_buf_close(&dIJAB);

  dpd_buf_init(&newtijab, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&dijab, CC_DENOM, 1, 6, 1, 6, 0, "dijab", 0, outfile);
  dpd_dirprd(&dijab, &newtijab, 0, outfile);
  dpd_buf_close(&newtijab);
  dpd_buf_close(&dijab);

  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_buf_init(&dIjAb, CC_DENOM, 0, 5, 0, 5, 0, "dIjAb", 0, outfile);
  dpd_dirprd(&dIjAb, &newtIjAb, 0, outfile);
  dpd_buf_close(&newtIjAb);
  dpd_buf_close(&dIjAb);
}
