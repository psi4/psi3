#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void WmbejL2(void)
{
  dpdbuf4 newL2, L2, W, Z, Z2;

  /* RHS += P(ij)P(ab)Limae * Wjebm */

  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,JB)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAJB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_sort(&Z, CC_TMP1, rqps, 10, 10, "Z(JA,IB)");
  dpd_buf4_sort(&Z, CC_TMP2, psrq, 10, 10, "Z(IB,JA)");
  dpd_buf4_sort(&Z, CC_TMP3, rspq, 10, 10, "Z(JB,IA)");
  dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(JA,IB)");
  dpd_buf4_axpy(&Z2, &Z, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 10, 10, 10, 0, "Z(IB,JA)");
  dpd_buf4_axpy(&Z2, &Z, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(JB,IA)");
  dpd_buf4_axpy(&Z2, &Z, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(IJ,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(IJ,AB)");
  dpd_buf4_init(&newL2, CC_LAMPS, 0, 0, 5, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&Z, &newL2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newL2);

  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ia,jb)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "Liajb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LiaJB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_sort(&Z, CC_TMP1, rqps, 10, 10, "Z(ja,ib)");
  dpd_buf4_sort(&Z, CC_TMP2, psrq, 10, 10, "Z(ib,ja)");
  dpd_buf4_sort(&Z, CC_TMP3, rspq, 10, 10, "Z(jb,ia)");
  dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ja,ib)");
  dpd_buf4_axpy(&Z2, &Z, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 10, 10, 10, 0, "Z(ib,ja)");
  dpd_buf4_axpy(&Z2, &Z, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(jb,ia)");
  dpd_buf4_axpy(&Z2, &Z, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(ij,ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(ij,ab)");
  dpd_buf4_init(&newL2, CC_LAMPS, 0, 0, 5, 2, 7, 0, "New Lijab");
  dpd_buf4_axpy(&Z, &newL2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newL2);


  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAJB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&W);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
  dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "Liajb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_contract444(&W, &L2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&W);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LiaJB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_contract444(&W, &L2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_init(&newL2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_buf4_axpy(&Z, &newL2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newL2);

  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ib,jA)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIbjA");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_contract444(&W, &L2, &Z, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&W);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LjAIb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
  dpd_contract444(&L2, &W, &Z, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&W);
  dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(Ij,bA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,bA)");
  dpd_buf4_sort(&Z, CC_TMP0, pqsr, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_init(&newL2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_buf4_axpy(&Z, &newL2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newL2);
}
