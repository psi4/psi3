#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void V_build(void)
{
  struct dpdbuf V, L, T;

  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "VMNIJ", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "Vmnij", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 0, 0, 0, 0, 0, "VMnIj", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAJB", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAJB", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "Viajb", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "Liajb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_contract222(&T, &L, &V, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAjb", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "Liajb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViaJB", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAJB", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_close(&T);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViAjB", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LIbjA", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_close(&L);
  dpd_buf_close(&V);

  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIaJb", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 10, 10, 10, 10, 0, "LjAIb", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_contract222(&T, &L, &V, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_close(&L);
  dpd_buf_close(&V);
}
