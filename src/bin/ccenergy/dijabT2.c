#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void dijabT2(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 dIJAB, dijab, dIjAb;

  dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
  dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
  dpd_buf4_dirprd(&dIJAB, &newtIJAB);
  dpd_buf4_close(&newtIJAB);
  dpd_buf4_close(&dIJAB);

  dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
  dpd_buf4_init(&dijab, CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
  dpd_buf4_dirprd(&dijab, &newtijab);
  dpd_buf4_close(&newtijab);
  dpd_buf4_close(&dijab);

  dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&dIjAb, &newtIjAb);
  dpd_buf4_close(&newtIjAb);
  dpd_buf4_close(&dIjAb);
}
