#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Gabcd(void)
{
  struct dpdbuf G, L, T;

  dpd_buf_init(&G, CC_GAMMA, 7, 7, 7, 7, 0, "GABCD", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_contract222(&L, &T, &G, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_close(&L);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 7, 7, 7, 7, 0, "Gabcd", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_contract222(&L, &T, &G, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_close(&L);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_contract222(&L, &T, &G, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_close(&L);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
}
