#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FmiL2(void)
{
  struct dpdbuf Lijab, LIJAB, LIjAb;
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct oe_dpdfile LFmit2, LFMIt2;

  /* RHS -= P(ij)*Limab*Fjm */
  dpd_oe_file_init(&LFMIt2, CC_OEI, 0, 0, "FMIt", 0, outfile);
  dpd_oe_file_init(&LFmit2, CC_OEI, 0, 0, "Fmit", 0, outfile);

  dpd_buf_init(&LIJAB, CC_LAMPS, 0, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&newLIJAB, CC_LAMPS, 0, 7, 2, 7, 0, "New LIJAB",
              0, outfile);
  dpd_contract221(&LIJAB, &LFMIt2, &newLIJAB, 1, 1, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&LFMIt2, &LIJAB, &newLIJAB, 1, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIJAB);
  dpd_buf_close(&newLIJAB);

  dpd_buf_init(&Lijab, CC_LAMPS, 0, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&newLijab, CC_LAMPS, 0, 7, 2, 7, 0, "New Lijab",
              0, outfile);
  dpd_contract221(&Lijab, &LFmit2, &newLijab, 1, 1, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&LFmit2, &Lijab, &newLijab, 1, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Lijab);
  dpd_buf_close(&newLijab);

  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
              0, outfile);
  dpd_contract221(&LIjAb, &LFmit2, &newLIjAb, 1, 1, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&LFMIt2, &LIjAb, &newLIjAb, 1, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIjAb);
  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&LFmit2);
  dpd_oe_file_close(&LFMIt2);
}
