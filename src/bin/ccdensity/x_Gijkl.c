#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* Gijkl(): computes non-R0 parts of Gijkl 2pdm */
/* Dijkl = 0.5 Lklef * Rijef + P(ij) Lklef * Rie * tjf */
/* Dijkl = R2L2_OOOO + P(ij) Lklef * Rie * tjf */

void x_Gijkl(void)
{
  dpdfile2 R1, T1;
  dpdbuf4 L2, I2, GIJKL, Gijkl, GIjKl;
  int L_irr, R_irr, G_irr;
  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  /* Gijkl += R2L2_OOOO */
  dpd_buf4_init(&GIJKL, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
  dpd_buf4_axpy(&I2, &GIJKL, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&GIJKL);

  dpd_buf4_init(&Gijkl, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_oooo");
  dpd_buf4_axpy(&I2, &Gijkl, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&Gijkl);

  dpd_buf4_init(&GIjKl, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_buf4_axpy(&I2, &GIjKl, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&GIjKl);

  /* GIJKL += - (Lklfe rie) tjf */
  dpd_buf4_init(&GIJKL, CC_GAMMA, G_irr, 0, 2, 2, 2, 0, "GIJKL");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I2, &T1, &GIJKL, 3, 1, 1, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &I2, &GIJKL, 1, 2, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&GIJKL);

  dpd_buf4_init(&Gijkl, CC_GAMMA, G_irr, 0, 2, 2, 2, 0, "Gijkl");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I2, &T1, &Gijkl, 3, 1, 1, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &I2, &Gijkl, 1, 2, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&Gijkl);

  dpd_buf4_init(&GIjKl, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I2, &T1, &GIjKl, 3, 1, 1, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &I2, &GIjKl, 1, 2, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&GIjKl);
}
