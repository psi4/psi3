#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void dijabL2(void)
{
  struct dpdbuf newLIJAB, newLijab, newLIjAb;
  struct dpdbuf dIJAB, dijab, dIjAb;

  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "New LIJAB", 0, outfile);
  dpd_buf_init(&dIJAB, CC_DENOM, 1, 6, 1, 6, 0, "dIJAB", 0, outfile);
  dpd_dirprd(&dIJAB, &newLIJAB, 0, outfile);
  dpd_buf_close(&newLIJAB);
  dpd_buf_close(&dIJAB);

  dpd_buf_init(&newLijab, CC_LAMPS, 2, 7, 2, 7, 0, "New Lijab", 0, outfile);
  dpd_buf_init(&dijab, CC_DENOM, 1, 6, 1, 6, 0, "dijab", 0, outfile);
  dpd_dirprd(&dijab, &newLijab, 0, outfile);
  dpd_buf_close(&newLijab);
  dpd_buf_close(&dijab);

  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb", 0, outfile);
  dpd_buf_init(&dIjAb, CC_DENOM, 0, 5, 0, 5, 0, "dIjAb", 0, outfile);
  dpd_dirprd(&dIjAb, &newLIjAb, 0, outfile);
  dpd_buf_close(&newLIjAb);
  dpd_buf_close(&dIjAb);
}

