#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double rohf_energy(void)
{
  double tIA_energy, tia_energy, tauIJAB_energy, tauijab_energy, tauIjAb_energy;
  dpdfile2 tIA, tia, fIA, fia;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, D;

  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
/*  dpd_file2_print(&tIA, outfile);  */
  tIA_energy = dpd_file2_dot(&fIA, &tIA);
  dpd_file2_close(&fIA);
  dpd_file2_close(&tIA);

  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
/*  dpd_file2_print(&tia, outfile); */
  tia_energy = dpd_file2_dot(&fia, &tia);
  dpd_file2_close(&fia);
  dpd_file2_close(&tia);

  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
/*  dpd_buf4_print(&tauIJAB, outfile);  */
  tauIJAB_energy = dpd_buf4_dot(&D, &tauIJAB);
  dpd_buf4_close(&tauIJAB);
  dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
/*  dpd_buf4_print(&tauijab, outfile); */
  tauijab_energy = dpd_buf4_dot(&D, &tauijab);
  dpd_buf4_close(&tauijab);
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
/*  dpd_buf4_print(&tauIjAb, outfile);  */
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

  return (tIA_energy + tia_energy +
	  tauIJAB_energy + tauijab_energy + tauIjAb_energy);
}
