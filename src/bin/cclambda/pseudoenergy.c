#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double pseudoenergy(void)
{
  double LIJAB_energy, Lijab_energy, LIjAb_energy;
  double LIA_energy=0.0, Lia_energy=0.0;
  struct dpdbuf LIJAB, Lijab, LIjAb, D;
  struct oe_dpdfile Lia, LIA, Fme, FME;

/*
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_oe_file_init(&FME, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_oe_file_init(&Lia, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&LIA, CC_OEI, 0, 1, "LIA", 0, outfile);

  LIA_energy = dpd_oe_dot(&FME,&LIA,0,outfile);
  Lia_energy = dpd_oe_dot(&Fme,&Lia,0,outfile);

  dpd_oe_file_close(&Lia);
  dpd_oe_file_close(&LIA);
  dpd_oe_file_close(&Fme);
  dpd_oe_file_close(&FME);
*/

  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
	       0, outfile);
  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  LIJAB_energy = dpd_dot(&D, &LIJAB, 0, outfile);
  dpd_buf_close(&LIJAB);
  dpd_buf_init(&Lijab, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  Lijab_energy = dpd_dot(&D, &Lijab, 0, outfile);
  dpd_buf_close(&Lijab);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  LIjAb_energy = dpd_dot(&D, &LIjAb, 0, outfile);
  dpd_buf_close(&LIjAb);
  dpd_buf_close(&D);

/*
  fprintf(outfile, "One A Energy = %20.14f\n", LIA_energy);
  fprintf(outfile, "One B Energy = %20.14f\n", Lia_energy);
*/

/*
  fprintf(outfile, "Two AA Energy = %20.14f\n", LIJAB_energy);
  fprintf(outfile, "Two BB Energy = %20.14f\n", Lijab_energy);
  fprintf(outfile, "Two AB Energy = %20.14f\n", LIjAb_energy);
*/

  return (LIA_energy+Lia_energy+LIJAB_energy + Lijab_energy + LIjAb_energy);
}
