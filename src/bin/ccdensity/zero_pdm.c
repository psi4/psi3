#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void zero_onepdm(void)
{
  dpdfile2 D;
  int G_irr;
  G_irr = params.G_irr;

  dpd_file2_init(&D, CC_OEI, G_irr, 0, 0, "DIJ");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, G_irr, 0, 0, "Dij");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, G_irr, 1, 1, "DAB");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, G_irr, 1, 1, "Dab");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, "DIA");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, "Dia");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, "DAI");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, "Dai");
  dpd_file2_scm(&D, 0.0);
  dpd_file2_close(&D);
}

void zero_twopdm(void) 
{
  dpdbuf4 G;
  int G_irr;
  G_irr = params.G_irr;

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);



  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_scm(&G, 0.0);
  dpd_buf4_close(&G);
}
