#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WijmbL2(void)
{
  dpdfile2 LIA, Lia;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
  dpdbuf4 X1, X2;

  /* RHS += -P(ab) Lma * Wijmb */
  dpd_file2_init(&LIA, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_OEI, 0, 0, 1, "Lia");

  dpd_buf4_init(&WMNIE, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE");
  dpd_buf4_init(&X1, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 1");
  dpd_contract424(&WMNIE, &LIA, &X1, 3, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&WMNIE);
  dpd_buf4_sort(&X1, CC_TMP1, pqsr, 2, 5, "X(2,5) 2");
  dpd_buf4_init(&X2, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 2");
  dpd_buf4_axpy(&X2, &X1, -1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&X1, &newLIJAB, 1.0);
  dpd_buf4_close(&newLIJAB);

  /*
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");

  dpd_buf4_init(&WMNIE, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE");
  dpd_buf4_sort(&WMNIE, CC_TMP0, pqsr, 2, 10, "WMNIE (M>N,IE)");
  dpd_contract424(&WMNIE, &LIA, &newLIJAB, 3, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WMNIE);

  dpd_buf4_init(&WMNIE, CC_TMP0, 0, 2, 10, 2, 10, 0, "WMNIE (M>N,IE)");
  dpd_contract244(&LIA, &WMNIE, &newLIJAB, 0, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&WMNIE);

  dpd_buf4_close(&newLIJAB);
  */

  dpd_buf4_init(&Wmnie, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie");
  dpd_buf4_init(&X1, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 1");
  dpd_contract424(&Wmnie, &Lia, &X1, 3, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&Wmnie);
  dpd_buf4_sort(&X1, CC_TMP1, pqsr, 2, 5, "X(2,5) 2");
  dpd_buf4_init(&X2, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 2");
  dpd_buf4_axpy(&X2, &X1, -1.0);
  dpd_buf4_close(&X2);
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New Lijab");
  dpd_buf4_axpy(&X1, &newLijab, 1.0);
  dpd_buf4_close(&newLijab);

  /*
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New Lijab");

  dpd_buf4_init(&Wmnie, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie");
  dpd_buf4_sort(&Wmnie, CC_TMP0, pqsr, 2, 10, "Wmnie (m>n,ie)");
  dpd_contract424(&Wmnie, &Lia, &newLijab, 3, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Wmnie);

  dpd_buf4_init(&Wmnie, CC_TMP0, 0, 2, 10, 2, 10, 0, "Wmnie (m>n,ie)");
  dpd_contract244(&Lia, &Wmnie, &newLijab, 0, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&Wmnie);

  dpd_buf4_close(&newLijab);
  */

  dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");

  dpd_buf4_init(&WMnIe, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_buf4_sort(&WMnIe, CC_TMP0, pqsr, 0, 10, "WMnIe (Mn,Ie)");
  dpd_buf4_close(&WMnIe);

  dpd_buf4_init(&WMnIe, CC_TMP0, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
  dpd_contract244(&LIA, &WMnIe, &newLIjAb, 0, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&WMnIe);

  dpd_buf4_init(&WmNiE, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE");
  dpd_buf4_sort(&WmNiE, CC_TMP0, qprs, 0, 11, "WmNiE (Nm,Ei)");
  dpd_buf4_close(&WmNiE);

  /* W(Nm,Ei) * L(i,b) --> L(Nm,Eb) */
  dpd_buf4_init(&WmNiE, CC_TMP0, 0, 0, 11, 0, 11, 0, "WmNiE (Nm,Ei)");
  dpd_contract424(&WmNiE, &Lia, &newLIjAb, 3, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WmNiE);

  dpd_buf4_close(&newLIjAb);

  dpd_file2_close(&Lia);
  dpd_file2_close(&LIA);
}

