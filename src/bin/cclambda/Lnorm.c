#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

extern double pseudoenergy(int L_irr);

void Lnorm(int L_irr, int R_index)
{
  dpdfile2 R1, L1, LIA, Lia, RIA, Ria;
  dpdbuf4 R2, L2, LIJAB, Lijab, LIjAb, RIJAB, Rijab, RIjAb;
  double tval, overlap, overlap0, overlap1, overlap2;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];

  sprintf(R1A_lbl, "RIA %d %d", L_irr, R_index);
  sprintf(R1B_lbl, "Ria %d %d", L_irr, R_index);
  sprintf(R2AA_lbl, "RIJAB %d %d", L_irr, R_index);
  sprintf(R2BB_lbl, "Rijab %d %d", L_irr, R_index);
  sprintf(R2AB_lbl, "RIjAb %d %d", L_irr, R_index);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    overlap0 = params.L0 * params.R0[L_irr][R_index];
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");

    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    overlap1 = dpd_file2_dot(&LIA, &R1);
    dpd_file2_close(&R1);
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
    overlap1 += dpd_file2_dot(&Lia, &R1);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    overlap2 = dpd_buf4_dot(&LIJAB, &R2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
    overlap2 += dpd_buf4_dot(&Lijab, &R2);
    dpd_buf4_close(&R2);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    overlap2 += dpd_buf4_dot(&LIjAb, &R2);
    dpd_buf4_close(&R2);

    overlap = overlap0 + overlap1 + overlap2;

    fprintf(outfile,"\n\tInitial  <L|R>  =     %15.10lf\n", overlap);

    dpd_file2_scm(&LIA, 1.0/overlap);
    dpd_file2_scm(&Lia, 1.0/overlap);
    dpd_buf4_scm(&LIJAB, 1.0/overlap);
    dpd_buf4_scm(&Lijab, 1.0/overlap);
    dpd_buf4_scm(&LIjAb, 1.0/overlap);

    fprintf(outfile,"\tNormalizing L...\n");
    fprintf(outfile,"\tL0 * R0 =     %15.10lf\n", overlap0/overlap);
    fprintf(outfile,"\tL1 * R1 =     %15.10lf\n", overlap1/overlap);
    fprintf(outfile,"\tL2 * R2 =     %15.10lf\n", overlap2/overlap);
    fprintf(outfile,"\t <L|R>  =     %15.10lf\n", overlap/overlap);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }

  tval = pseudoenergy(L_irr);
  fprintf(outfile,"\tPseudoenergy w/normed L = %15.10lf\n",tval);

  return;
}
