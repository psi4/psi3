#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

double rhf_energy(void)
{
  double tIA_energy, tauIjAb_energy;
  dpdfile2 tIA, fIA;
  dpdbuf4 tauIjAb, D;

  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  tIA_energy = 2 * dpd_file2_dot(&fIA, &tIA);
  dpd_file2_close(&fIA);
  dpd_file2_close(&tIA);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ji|ab>");
  dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  tauIjAb_energy = dpd_buf4_dot(&D, &tauIjAb);
  dpd_buf4_close(&tauIjAb);
  dpd_buf4_close(&D);

  /*
    fprintf(outfile, "One A Energy = %20.14f\n", tIA_energy);
    fprintf(outfile, "One B Energy = %20.14f\n", tia_energy);
    fprintf(outfile, "Two AA Energy = %20.14f\n", tauIJAB_energy);
    fprintf(outfile, "Two BB Energy = %20.14f\n", tauijab_energy);
    fprintf(outfile, "Two AB Energy = %20.14f\n", tauIjAb_energy);
  */

  return (tIA_energy + tauIjAb_energy);
}
