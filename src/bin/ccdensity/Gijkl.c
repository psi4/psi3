#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Gijkl(void)
{
  dpdbuf4 V, G;
  int G_irr;
  G_irr = params.G_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&V, CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
    dpd_buf4_copy(&V, CC_GAMMA, "GIJKL");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }

    dpd_buf4_init(&V, CC_MISC, G_irr, 2, 2, 2, 2, 0, "Vmnij");
    dpd_buf4_copy(&V, CC_GAMMA, "Gijkl");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }

    dpd_buf4_init(&V, CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
    dpd_buf4_copy(&V, CC_GAMMA, "GIjKl");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&V, CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
    dpd_buf4_copy(&V, CC_GAMMA, "GIJKL");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }

    dpd_buf4_init(&V, CC_MISC, G_irr, 12, 12, 12, 12, 0, "Vmnij");
    dpd_buf4_copy(&V, CC_GAMMA, "Gijkl");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 12, 12, 12, 12, 0, "Gijkl");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }

    dpd_buf4_init(&V, CC_MISC, G_irr, 22, 22, 22, 22, 0, "VMnIj");
    dpd_buf4_copy(&V, CC_GAMMA, "GIjKl");
    dpd_buf4_close(&V);
    if (params.ground) {
      dpd_buf4_init(&G, CC_GAMMA, G_irr, 22, 22, 22, 22, 0, "GIjKl");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
    }
  }
}
