#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

double rhf_energy(void);
double rohf_energy(void);
double uhf_energy(void);

double energy(void)
{
  if(params.ref == 0) return(rhf_energy());
  else if(params.ref == 1) return(rohf_energy());
  else if(params.ref == 2) return(uhf_energy());
}

double rhf_energy(void)
{
  double tauIjAb_energy;
  dpdbuf4 tauIjAb, D;

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  tauIjAb_energy = dpd_buf4_dot(&D, &tauIjAb);
  dpd_buf4_close(&tauIjAb);
  dpd_buf4_close(&D);

  /*
    fprintf(outfile, "Two AB Energy = %20.14f\n", tauIjAb_energy);
  */

  return (tauIjAb_energy);
}

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

double uhf_energy(void)
{
  double E2AA, E2BB, E2AB;
  dpdbuf4 T2, D;

  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);

  return(E2AA + E2BB + E2AB);
}
