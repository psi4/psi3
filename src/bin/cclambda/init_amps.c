#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void init_amps(int L_irr, int R_index)
{
  double norm;
  dpdfile2 T1, R1, LIA, Lia;
  dpdbuf4 T2, R2, LIJAB, Lijab, LIjAb;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];


  /* Restart from previous amplitudes if we can/should */
  /* Need to adjust this for new I/O
     if(params.restart && flen(CC_LIA) && flen(CC_Lia) && flen(CC_LIJAB)
     && flen(CC_Lijab) && flen(CC_LIjAb)) return;
  */

  /* ground state guess L <= T */
  /* excited state guess L <= R0 * T + R */
  if (params.ground || L_irr == 0) {
    if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_copy(&T1, CC_OEI, "LIA");
      dpd_file2_close(&T1);
  
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
      dpd_file2_copy(&T1, CC_OEI, "Lia");
      dpd_file2_close(&T1);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_buf4_copy(&T2, CC_LAMBDA, "LIJAB");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
      dpd_buf4_copy(&T2, CC_LAMBDA, "Lijab");
      dpd_buf4_close(&T2);

      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_copy(&T2, CC_LAMBDA, "LIjAb");
      dpd_buf4_close(&T2);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_copy(&T1, CC_OEI, "LIA");
      dpd_file2_close(&T1);
  
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_copy(&T1, CC_OEI, "Lia");
      dpd_file2_close(&T1);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_buf4_copy(&T2, CC_LAMBDA, "LIJAB");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
      dpd_buf4_copy(&T2, CC_LAMBDA, "Lijab");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      dpd_buf4_copy(&T2, CC_LAMBDA, "LIjAb");
      dpd_buf4_close(&T2);
    }
  }

  if (!params.ground) {
    sprintf(R1A_lbl, "RIA %d %d", L_irr, R_index);
    sprintf(R1B_lbl, "Ria %d %d", L_irr, R_index);
    sprintf(R2AA_lbl, "RIJAB %d %d", L_irr, R_index);
    sprintf(R2BB_lbl, "Rijab %d %d", L_irr, R_index);
    sprintf(R2AB_lbl, "RIjAb %d %d", L_irr, R_index);

    /* multiply by R0 and create nonsymmetric L files */
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    if (params.ref <= 1) {
      dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    }
    else {
      dpd_file2_init(&Lia, CC_OEI, L_irr, 2, 3, "Lia");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    }

    dpd_file2_scm(&LIA, params.R0[L_irr][R_index]);
    dpd_file2_scm(&Lia, params.R0[L_irr][R_index]);
    dpd_buf4_scm(&LIJAB, params.R0[L_irr][R_index]);
    dpd_buf4_scm(&Lijab, params.R0[L_irr][R_index]);
    dpd_buf4_scm(&LIjAb, params.R0[L_irr][R_index]);
  
      /* add R1 and R2 */
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    dpd_file2_axpy(&R1, &LIA, 1.0, 0);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    dpd_buf4_axpy(&R2, &LIJAB, 1.0);
    dpd_buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      dpd_file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      dpd_buf4_axpy(&R2, &Lijab, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      dpd_buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      dpd_file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      dpd_buf4_axpy(&R2, &Lijab, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      dpd_buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_buf4_close(&R2);
    }
  
    /* dot L and R together */
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_file2_dot(&LIA, &R1);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_buf4_dot(&LIJAB, &R2);
    dpd_buf4_close(&R2);
    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
  
    fprintf(outfile,"\tInitial overlap of initial guess <R|L> = %15.10lf\n", norm);
  
    dpd_file2_scm(&LIA, 1.0/norm);
    dpd_file2_scm(&Lia, 1.0/norm);
    dpd_buf4_scm(&LIJAB, 1.0/norm);
    dpd_buf4_scm(&Lijab, 1.0/norm);
    dpd_buf4_scm(&LIjAb, 1.0/norm);
  
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_file2_dot(&LIA, &R1);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_buf4_dot(&LIJAB, &R2);
    dpd_buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    fprintf(outfile,"\tChecking overlap of initial guess <R|L> = %15.10lf\n", norm);
  
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }

#ifdef EOM_DEBUG
  fprintf(outfile,"initial guess\n");
  dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
  dpd_file2_print(&LIA,outfile);
  dpd_file2_close(&LIA);
#endif
}
