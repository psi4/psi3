#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double LR_overlap_ROHF(int IRR, int L_index, int R_index);
double LR_overlap_RHF(int IRR, int L_index, int R_index);

void check_ortho(void) {
  int L_state_index, root_L_irr, L_irr, nstates=0;
  int R_state_index, root_R_irr, R_irr;
  double **O, tval;

  for (L_irr=0; L_irr<moinfo.nirreps; ++L_irr)
    nstates += params.Ls_per_irrep[L_irr];

  if (params.ref <= 1) {
    L_state_index = 0;
    O = block_matrix(nstates,nstates);
    for (L_irr=0; L_irr<moinfo.nirreps; ++L_irr)
      for (root_L_irr=0; root_L_irr< params.Ls_per_irrep[L_irr]; ++root_L_irr) {
        ++L_state_index;
  
        R_state_index = 0;
        for (R_irr=0; R_irr<moinfo.nirreps; ++R_irr)
          for (root_R_irr=0; root_R_irr< params.Ls_per_irrep[R_irr]; ++root_R_irr) {
            ++R_state_index;
  
            if (L_irr == R_irr)
              tval = LR_overlap_ROHF(L_irr, root_L_irr, root_R_irr);
            else
              tval = -99.0;
  
            O[L_state_index-1][R_state_index-1] = tval;
          }
      }
  
    fprintf(outfile,"\t<L|R> overlap matrix with ROHF quantities (-99 => 0 by symmetry\n");
    print_mat(O, nstates, nstates, outfile);
    free_block(O);
  }

  if (params.ref == 0) { /* test RHF quantities */
    O = block_matrix(nstates, nstates);
    L_state_index = 0;
    for (L_irr=0; L_irr<moinfo.nirreps; ++L_irr)
      for (root_L_irr=0; root_L_irr< params.Ls_per_irrep[L_irr]; ++root_L_irr) {
        ++L_state_index;
 
        R_state_index = 0;
        for (R_irr=0; R_irr<moinfo.nirreps; ++R_irr)
          for (root_R_irr=0; root_R_irr< params.Ls_per_irrep[R_irr]; ++root_R_irr) {
            ++R_state_index;
 
            if (L_irr == R_irr)
              tval = LR_overlap_RHF(L_irr, root_L_irr, root_R_irr);
            else
              tval = -99.0;
 
            O[L_state_index-1][R_state_index-1] = tval;
          }
      }
 
    fprintf(outfile,"\t<L|R> overlap matrix with RHF quantities (-99 => 0 by symmetry)\n"); 
    print_mat(O, nstates, nstates, outfile);
    free_block(O);
  }
  return;
}

double LR_overlap_ROHF(int IRR, int L_index, int R_index) {
  double overlap;
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  char L1A_lbl[32], L1B_lbl[32], L2AA_lbl[32], L2BB_lbl[32], L2AB_lbl[32];

  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
  sprintf(R1B_lbl, "Ria %d %d", IRR, R_index);
  sprintf(R2AA_lbl, "RIJAB %d %d", IRR, R_index);
  sprintf(R2BB_lbl, "Rijab %d %d", IRR, R_index);
  sprintf(R2AB_lbl, "RIjAb %d %d", IRR, R_index);

  sprintf(L1A_lbl, "LIA %d %d", IRR, L_index);
  sprintf(L1B_lbl, "Lia %d %d", IRR, L_index);
  sprintf(L2AA_lbl, "LIJAB %d %d", IRR, L_index);
  sprintf(L2BB_lbl, "Lijab %d %d", IRR, L_index);
  sprintf(L2AB_lbl, "LIjAb %d %d", IRR, L_index);

  dpd_file2_init(&R1, CC_RAMPS, IRR, 0, 1, R1A_lbl);
  dpd_file2_init(&L1, CC_LAMPS, IRR, 0, 1, L1A_lbl);
  overlap = dpd_file2_dot(&L1, &R1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_file2_init(&R1, CC_RAMPS, IRR, 0, 1, R1B_lbl);
  dpd_file2_init(&L1, CC_LAMPS, IRR, 0, 1, L1B_lbl);
  overlap += dpd_file2_dot(&L1, &R1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_buf4_init(&R2, CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2AA_lbl);
  dpd_buf4_init(&L2, CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
  overlap += dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&R2, CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2BB_lbl);
  dpd_buf4_init(&L2, CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2BB_lbl);
  overlap += dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&R2, CC_RAMPS, IRR, 0, 5, 0, 5, 0, R2AB_lbl);
  dpd_buf4_init(&L2, CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
  overlap += dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);
  
  return overlap;
}

double LR_overlap_RHF(int IRR, int L_index, int R_index) {
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  double overlap, overlap2, overlap3;
  char L1A_lbl[32], R1A_lbl[32], lbl[32];

  sprintf(L1A_lbl, "LIA %d %d", IRR, L_index);
  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);

  dpd_file2_init(&R1, CC_RAMPS, IRR, 0, 1, R1A_lbl);
  dpd_file2_init(&L1, CC_LAMPS, IRR, 0, 1, L1A_lbl);
  overlap = 2.0 * dpd_file2_dot(&L1, &R1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  sprintf(lbl, "2RIjAb - RIjbA %d %d", IRR, R_index);
  dpd_buf4_init(&R2, CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "LIjAb %d %d", IRR, L_index);
  dpd_buf4_init(&L2, CC_LAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  overlap2 = dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&R2);

  sprintf(lbl, "2LIjAb - LIjbA %d %d", IRR, L_index);
  dpd_buf4_init(&L2, CC_LAMPS, IRR, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "RIjAb %d %d", IRR, R_index);
  dpd_buf4_init(&R2, CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  overlap3 = dpd_buf4_dot(&L2, &R2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  if (fabs(overlap2 - overlap3) > 1E-14) {
    fprintf(outfile,"Bad anti-symmetry detected in RHF quantities\n");
    fprintf(outfile,"error: %15.10lf\n",overlap2-overlap3);
  }

  overlap += overlap2;
  return overlap;
}
