#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WejabL2(void)
{
  struct dpdbuf Wamef, WAmEf, WaMeF, WAMEF;
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf Ltmp;
  struct oe_dpdfile LIA, Lia;
  
  /* RHS += P(ij) Lie * Wejab */
  dpd_oe_file_init(&LIA, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&Lia, CC_OEI, 0, 1, "Lia", 0, outfile);

  dpd_buf_init(&WAMEF, CC_HBAR, 10, 7, 10, 7, 0, "WAMEF", 0, outfile);
  dpd_buf_init(&newLIJAB, CC_LAMPS, 0, 7, 2, 7, 0, "New LIJAB",
              0, outfile);
  dpd_contract221(&WAMEF, &LIA, &newLIJAB, 1, 1, 1, -1.0, 1.0, 0, outfile);

  dpd_buf_init(&Ltmp, CC_TMP0, 0, 7, 0, 7, 0, "LIJAB (JI,A>B)", 0, outfile);
  dpd_contract221(&WAMEF, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_swap12(&Ltmp, CC_TMP1, 0, 7, "LIJAB (IJ,A>B)", 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_init(&Ltmp, CC_TMP1, 0, 7, 0, 7, 0, "LIJAB (IJ,A>B)", 0, outfile);
  dpd_axpy(&Ltmp, &newLIJAB, 1.0, 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&newLIJAB);
  dpd_buf_close(&WAMEF);

  dpd_buf_init(&Wamef, CC_HBAR, 10, 7, 10, 7, 0, "Wamef", 0, outfile);
  dpd_buf_init(&newLijab, CC_LAMPS, 0, 7, 2, 7, 0, "New Lijab",
              0, outfile);
  dpd_contract221(&Wamef, &Lia, &newLijab, 1, 1, 1, -1.0, 1.0, 0, outfile);

  dpd_buf_init(&Ltmp, CC_TMP0, 0, 7, 0, 7, 0, "Lijab (ji,a>b)", 0, outfile);
  dpd_contract221(&Wamef, &Lia, &Ltmp, 1, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_swap12(&Ltmp, CC_TMP1, 0, 7, "Lijab (ij,a>b)", 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_init(&Ltmp, CC_TMP1, 0, 7, 0, 7, 0, "Lijab (ij,a>b)", 0, outfile);
  dpd_axpy(&Ltmp, &newLijab, 1.0, 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&newLijab);
  dpd_buf_close(&Wamef);

  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
             0, outfile);

  dpd_buf_init(&WaMeF, CC_HBAR, 10, 5, 10, 5, 0, "WaMeF", 0, outfile);
  dpd_swap34(&WaMeF, CC_TMP0, 10, 5, "WaMeF (Ma,Fe)", 0, outfile);
  dpd_buf_close(&WaMeF);

  dpd_buf_init(&WaMeF, CC_TMP0, 10, 5, 10, 5, 0, "WaMeF (Ma,Fe)", 0, outfile);
  dpd_contract221(&WaMeF, &Lia, &newLIjAb, 1, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WaMeF);

  dpd_buf_init(&WAmEf, CC_HBAR, 10, 5, 10, 5, 0, "WAmEf", 0, outfile);
  dpd_buf_init(&Ltmp, CC_TMP0, 0, 5, 0, 5, 0, "Ltmp (ji,ab)", 0, outfile);
  dpd_contract221(&WAmEf, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_swap12(&Ltmp, CC_TMP1, 0, 5, "Lijab (ij,ab)", 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&WAmEf);

  dpd_buf_init(&Ltmp, CC_TMP1, 0, 5, 0, 5, 0, "Lijab (ij,ab)", 0, outfile);
  dpd_axpy(&Ltmp, &newLIjAb, 1.0, 0, outfile);
  dpd_buf_close(&Ltmp);

  dpd_buf_close(&newLIjAb);
  dpd_oe_file_close(&Lia);
  dpd_oe_file_close(&LIA);
}
