#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double pseudoenergy(void)
{
  double LIJAB_energy, Lijab_energy, LIjAb_energy;
  double LIA_energy=0.0, Lia_energy=0.0;
  dpdbuf4 LIJAB, Lijab, LIjAb, D;
  dpdfile2 Lia, LIA, Fme, FME;

  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Lia, CC_OEI, 0, 0, 1, "Lia");
  dpd_file2_init(&LIA, CC_OEI, 0, 0, 1, "LIA");

  LIA_energy = dpd_file2_dot(&FME,&LIA);
  Lia_energy = dpd_file2_dot(&Fme,&Lia);

  dpd_file2_close(&Lia);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Fme);
  dpd_file2_close(&FME);

  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
  LIJAB_energy = dpd_buf4_dot(&D, &LIJAB);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
  Lijab_energy = dpd_buf4_dot(&D, &Lijab);
  dpd_buf4_close(&Lijab);
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  LIjAb_energy = dpd_buf4_dot(&D, &LIjAb);
  dpd_buf4_close(&LIjAb);
  dpd_buf4_close(&D);

  /*
  fprintf(outfile, "One A Energy = %20.14f\n", LIA_energy);
  fprintf(outfile, "One B Energy = %20.14f\n", Lia_energy);

  fprintf(outfile, "Two AA Energy = %20.14f\n", LIJAB_energy);
  fprintf(outfile, "Two BB Energy = %20.14f\n", Lijab_energy);
  fprintf(outfile, "Two AB Energy = %20.14f\n", LIjAb_energy);
  */

  return (LIJAB_energy + Lijab_energy + LIjAb_energy);
}
