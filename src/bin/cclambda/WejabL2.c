#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void WejabL2(void)
{
  dpdbuf4 Wamef, WAmEf, WaMeF, WAMEF;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Ltmp;
  dpdfile2 LIA, Lia;
  dpdbuf4 X1, X2;
  
  /* RHS += P(ij) Lie * Wejab */
  dpd_file2_init(&LIA, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_OEI, 0, 0, 1, "Lia");

  dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
  dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
  dpd_contract424(&WAMEF, &LIA, &X1, 1, 1, 1, -1.0, 0.0);
  dpd_buf4_close(&WAMEF);
  dpd_buf4_sort(&X1, CC_TMP1, qprs, 0, 7, "X(0,7) 2");
  dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
  dpd_buf4_axpy(&X2, &X1, -1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&X1, &newLIJAB, 1.0);
  dpd_buf4_close(&X1);
  dpd_buf4_close(&newLIJAB);

  /*
  dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New LIJAB");
  dpd_contract424(&WAMEF, &LIA, &newLIJAB, 1, 1, 1, -1.0, 1.0);

  dpd_buf4_init(&Ltmp, CC_TMP0, 0, 0, 7, 0, 7, 0, "LIJAB (JI,A>B)");
  dpd_contract424(&WAMEF, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0);
  dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 7, "LIJAB (IJ,A>B)");
  dpd_buf4_close(&Ltmp);
  dpd_buf4_init(&Ltmp, CC_TMP1, 0, 0, 7, 0, 7, 0, "LIJAB (IJ,A>B)");
  dpd_buf4_axpy(&Ltmp, &newLIJAB, 1.0);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&newLIJAB);
  dpd_buf4_close(&WAMEF);
  */

  dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
  dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
  dpd_contract424(&Wamef, &Lia, &X1, 1, 1, 1, -1.0, 0.0);
  dpd_buf4_close(&Wamef);
  dpd_buf4_sort(&X1, CC_TMP1, qprs, 0, 7, "X(0,7) 2");
  dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
  dpd_buf4_axpy(&X2, &X1, -1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New Lijab");
  dpd_buf4_axpy(&X1, &newLijab, 1.0);
  dpd_buf4_close(&X1);
  dpd_buf4_close(&newLijab);

  /*
  dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New Lijab");
  dpd_contract424(&Wamef, &Lia, &newLijab, 1, 1, 1, -1.0, 1.0);

  dpd_buf4_init(&Ltmp, CC_TMP0, 0, 0, 7, 0, 7, 0, "Lijab (ji,a>b)");
  dpd_contract424(&Wamef, &Lia, &Ltmp, 1, 1, 1, 1.0, 0.0);
  dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 7, "Lijab (ij,a>b)");
  dpd_buf4_close(&Ltmp);
  dpd_buf4_init(&Ltmp, CC_TMP1, 0, 0, 7, 0, 7, 0, "Lijab (ij,a>b)");
  dpd_buf4_axpy(&Ltmp, &newLijab, 1.0);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&newLijab);
  dpd_buf4_close(&Wamef);
  */

  dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");

  dpd_buf4_init(&WaMeF, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
  dpd_buf4_sort(&WaMeF, CC_TMP0, pqsr, 10, 5, "WaMeF (Ma,Fe)");
  dpd_buf4_close(&WaMeF);

  dpd_buf4_init(&WaMeF, CC_TMP0, 0, 10, 5, 10, 5, 0, "WaMeF (Ma,Fe)");
  dpd_contract424(&WaMeF, &Lia, &newLIjAb, 1, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&WaMeF);

  dpd_buf4_init(&WAmEf, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_buf4_init(&Ltmp, CC_TMP0, 0, 0, 5, 0, 5, 0, "Ltmp (ji,ab)");
  dpd_contract424(&WAmEf, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0);
  dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 5, "Lijab (ij,ab)");
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&WAmEf);

  dpd_buf4_init(&Ltmp, CC_TMP1, 0, 0, 5, 0, 5, 0, "Lijab (ij,ab)");
  dpd_buf4_axpy(&Ltmp, &newLIjAb, 1.0);
  dpd_buf4_close(&Ltmp);

  dpd_buf4_close(&newLIjAb);
  dpd_file2_close(&Lia);
  dpd_file2_close(&LIA);
}
