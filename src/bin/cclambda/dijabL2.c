#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void dijabL2(void)
{
  dpdbuf4 L2, newLIJAB, newLijab, newLIjAb;
  dpdbuf4 d2, dIJAB, dijab, dIjAb;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&dIJAB, &newLIJAB);
    dpd_buf4_close(&newLIJAB);
    dpd_buf4_close(&dIJAB);

    dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&dijab, CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
    dpd_buf4_dirprd(&dijab, &newLijab);
    dpd_buf4_close(&newLijab);
    dpd_buf4_close(&dijab);

    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&dIjAb, &newLIjAb);
    dpd_buf4_close(&newLIjAb);
    dpd_buf4_close(&dIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&d2, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&d2, CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&d2, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

  }
}

