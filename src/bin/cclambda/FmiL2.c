#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void FmiL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdfile2 LFmit2, LFMIt2;
  dpdbuf4 X1, X2;

  /* RHS -= P(ij)*Limab*Fjm */
  dpd_file2_init(&LFMIt2, CC_OEI, 0, 0, 0, "FMIt");
  dpd_file2_init(&LFmit2, CC_OEI, 0, 0, 0, "Fmit");

  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
  dpd_contract424(&LIJAB, &LFMIt2, &X1, 1, 1, 1, -1.0, 0.0);
  dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
  dpd_contract244(&LFMIt2, &LIJAB, &X2, 1, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_axpy(&X1, &X2, 1.0);
  dpd_buf4_close(&X1);
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_close(&newLIJAB);

  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
  dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
  dpd_contract424(&Lijab, &LFmit2, &X1, 1, 1, 1, -1.0, 0.0);
  dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
  dpd_contract244(&LFmit2, &Lijab, &X2, 1, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&Lijab);
  dpd_buf4_axpy(&X1, &X2, 1.0);
  dpd_buf4_close(&X1);
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New Lijab");
  dpd_buf4_axpy(&X2, &newLijab, 1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_close(&newLijab);

  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_contract424(&LIjAb, &LFmit2, &newLIjAb, 1, 1, 1, -1.0, 1.0);
  dpd_contract244(&LFMIt2, &LIjAb, &newLIjAb, 1, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&LIjAb);
  dpd_buf4_close(&newLIjAb);

  dpd_file2_close(&LFmit2);
  dpd_file2_close(&LFMIt2);
}
