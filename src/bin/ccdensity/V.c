#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void V_build(void)
{
  dpdbuf4 V, L, T;

  dpd_buf4_init(&V, CC_MISC, 0, 2, 2, 2, 2, 0, "VMNIJ");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&L, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 2, 2, 2, 2, 0, "Vmnij");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_init(&L, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 0, 0, 0, 0, 0, "VMnIj");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIAJB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAJB");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "Viajb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "Liajb");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIAjb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "Liajb");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAJB");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&T);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIbjA");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
  dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&L);
  dpd_buf4_close(&V);

  dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIaJb");
  dpd_buf4_init(&L, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LjAIb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&L);
  dpd_buf4_close(&V);
}
