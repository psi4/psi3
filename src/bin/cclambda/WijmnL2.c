#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WijmnL2(void)
{
  struct dpdbuf Lijab, LIJAB, LIjAb;
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf WMNIJ, Wmnij, WMnIj;

   /* RHS += Lmnab*Wijmn */
  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "New LIJAB",
              0, outfile);
  dpd_buf_init(&WMNIJ, CC_HBAR, 2, 2, 2, 2, 0, "WMNIJ", 0, outfile);
  dpd_contract222(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMNIJ);
  dpd_buf_close(&LIJAB);
  dpd_buf_close(&newLIJAB);

  dpd_buf_init(&Lijab, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&newLijab, CC_LAMPS, 2, 7, 2, 7, 0, "New Lijab",
              0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 2, 2, 2, 2, 0, "Wmnij", 0, outfile);
  dpd_contract222(&Wmnij, &Lijab, &newLijab, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wmnij);
  dpd_buf_close(&Lijab);
  dpd_buf_close(&newLijab);

  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
             0, outfile);
  dpd_buf_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);
  dpd_contract222(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMnIj);
  dpd_buf_close(&LIjAb);
  dpd_buf_close(&newLIjAb);
}

