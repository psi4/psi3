#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void x_Gijab_ROHF(void)
{
  int h, nirreps, i, a, m, e, I, A, M, E, Isym, Asym, Msym, Esym, row, col;
  int R_irr, L_irr, G_irr;
  double value;
  dpdfile2 T1, L1, g, ZZ, ZZ2, T1A, T1B;
  dpdbuf4 G, L, T, V, Z, Z1, Z2;

  nirreps = moinfo.nirreps;
  R_irr = params.R_irr; L_irr = params.L_irr; G_irr = params.G_irr;

  /* (1-R0) * Tau(IJ,AB) */
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_scmcopy(&T, EOM_TMP0, "GIJAB", 1.0 - params.R0);
  dpd_buf4_close(&T);

  /* (1-R0) * Tau(ij,ab) */
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_scmcopy(&T, EOM_TMP0, "Gijab", 1.0 - params.R0);
  dpd_buf4_close(&T);

  /* (1-R0) * Tau(Ij,Ab) */
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_scmcopy(&T, EOM_TMP0, "GIjAb", 1.0 - params.R0);
  dpd_buf4_close(&T);

  /* add to ground state parts */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  return;
}
