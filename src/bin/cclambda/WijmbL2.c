#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WijmbL2(void)
{
  struct oe_dpdfile LIA, Lia;
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf WMNIE, Wmnie, WMnIe, WmNiE;

  /* RHS += -P(ab) Lma * Wijmb */
  dpd_oe_file_init(&LIA, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&Lia, CC_OEI, 0, 1, "Lia", 0, outfile);

  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "New LIJAB",
              0, outfile);

  dpd_buf_init(&WMNIE, CC_HBAR, 2, 11, 2, 11, 0, "WMNIE", 0, outfile);
  dpd_swap34(&WMNIE, CC_TMP0, 2, 10, "WMNIE (M>N,IE)", 0, outfile);
  dpd_contract221(&WMNIE, &LIA, &newLIJAB, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMNIE);

  dpd_buf_init(&WMNIE, CC_TMP0, 2, 10, 2, 10, 0, "WMNIE (M>N,IE)", 0, outfile);
  dpd_contract212(&LIA, &WMNIE, &newLIJAB, 0, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMNIE);

  dpd_buf_close(&newLIJAB);


  dpd_buf_init(&newLijab, CC_LAMPS, 2, 5, 2, 7, 0, "New Lijab",
              0, outfile);

  dpd_buf_init(&Wmnie, CC_HBAR, 2, 11, 2, 11, 0, "Wmnie", 0, outfile);
  dpd_swap34(&Wmnie, CC_TMP0, 2, 10, "Wmnie (m>n,ie)", 0, outfile);
  dpd_contract221(&Wmnie, &Lia, &newLijab, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wmnie);

  dpd_buf_init(&Wmnie, CC_TMP0, 2, 10, 2, 10, 0, "Wmnie (m>n,ie)", 
               0, outfile);
  dpd_contract212(&Lia, &Wmnie, &newLijab, 0, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wmnie);

  dpd_buf_close(&newLijab);

  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
             0, outfile);

  dpd_buf_init(&WMnIe, CC_HBAR, 0, 11, 0, 11, 0, "WMnIe", 0, outfile);
  dpd_swap34(&WMnIe, CC_TMP0, 0, 10, "WMnIe (Mn,Ie)", 0, outfile);
  dpd_buf_close(&WMnIe);

  dpd_buf_init(&WMnIe, CC_TMP0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)", 0, outfile);
  dpd_contract212(&LIA, &WMnIe, &newLIjAb, 0, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMnIe);

  dpd_buf_init(&WmNiE, CC_HBAR, 0, 11, 0, 11, 0, "WmNiE", 0, outfile);
  dpd_swap12(&WmNiE, CC_TMP0, 0, 11, "WmNiE (Nm,Ei)", 0, outfile);
  dpd_buf_close(&WmNiE);

  /* W(Nm,Ei) * L(i,b) --> L(Nm,Eb) */
  dpd_buf_init(&WmNiE, CC_TMP0, 0, 11, 0, 11, 0, "WmNiE (Nm,Ei)", 0, outfile);
  dpd_contract221(&WmNiE, &Lia, &newLIjAb, 3, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WmNiE);

  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&Lia);
  dpd_oe_file_close(&LIA);
}

