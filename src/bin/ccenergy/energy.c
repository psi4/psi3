#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double energy(void)
{
  double tIA_energy, tia_energy, tauIJAB_energy, tauijab_energy, tauIjAb_energy;
  struct oe_dpdfile tIA, tia, fIA, fia;
  struct dpdbuf tauIJAB, tauijab, tauIjAb, D;

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  tIA_energy = dpd_oe_dot(&fIA, &tIA, 0, outfile);
  dpd_oe_file_close(&fIA);
  dpd_oe_file_close(&tIA);

  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  tia_energy = dpd_oe_dot(&fia, &tia, 0, outfile);
  dpd_oe_file_close(&fia);
  dpd_oe_file_close(&tia);

  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
	       0, outfile);
  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  tauIJAB_energy = dpd_dot(&D, &tauIJAB, 0, outfile);
  dpd_buf_close(&tauIJAB);
  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  tauijab_energy = dpd_dot(&D, &tauijab, 0, outfile);
  dpd_buf_close(&tauijab);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  tauIjAb_energy = dpd_dot(&D, &tauIjAb, 0, outfile);
  dpd_buf_close(&tauIjAb);
  dpd_buf_close(&D);

  /*
  fprintf(outfile, "One A Energy = %20.14f\n", tIA_energy);
  fprintf(outfile, "One B Energy = %20.14f\n", tia_energy);
  fprintf(outfile, "Two AA Energy = %20.14f\n", tauIJAB_energy);
  fprintf(outfile, "Two BB Energy = %20.14f\n", tauijab_energy);
  fprintf(outfile, "Two AB Energy = %20.14f\n", tauIjAb_energy);
  */

  return (tIA_energy + tia_energy +
	  tauIJAB_energy + tauijab_energy + tauIjAb_energy);
}
