#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Gijkl(void)
{
  struct dpdbuf V, G;

  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "VMNIJ", 0, outfile);
  dpd_copy(&V, CC_GAMMA, "GIJKL", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&G, CC_GAMMA, 2, 2, 2, 2, 0, "GIJKL", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);

  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "Vmnij", 0, outfile);
  dpd_copy(&V, CC_GAMMA, "Gijkl", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&G, CC_GAMMA, 2, 2, 2, 2, 0, "Gijkl", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);

  dpd_buf_init(&V, CC_MISC, 0, 0, 0, 0, 0, "VMnIj", 0, outfile);
  dpd_copy(&V, CC_GAMMA, "GIjKl", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, "GIjKl", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
}
