#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Gijkl(void)
{
  dpdbuf4 V, G;

  dpd_buf4_init(&V, CC_MISC, 0, 2, 2, 2, 2, 0, "VMNIJ");
  dpd_buf4_copy(&V, CC_GAMMA, "GIJKL");
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 2, 2, 2, 0, "GIJKL");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);

  dpd_buf4_init(&V, CC_MISC, 0, 2, 2, 2, 2, 0, "Vmnij");
  dpd_buf4_copy(&V, CC_GAMMA, "Gijkl");
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 2, 2, 2, 0, "Gijkl");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);

  dpd_buf4_init(&V, CC_MISC, 0, 0, 0, 0, 0, 0, "VMnIj");
  dpd_buf4_copy(&V, CC_GAMMA, "GIjKl");
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);
}
