#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

extern void Lnorm(void);

void init_amps(void)
{
  double norm;
  dpdfile2 T1, R1, LIA, Lia;
  dpdbuf4 T2, R2, LIJAB, Lijab, LIjAb;

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
      dpd_buf4_copy(&T2, CC_LAMPS, "LIJAB");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
      dpd_buf4_copy(&T2, CC_LAMPS, "Lijab");
      dpd_buf4_close(&T2);

      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_copy(&T2, CC_LAMPS, "LIjAb");
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
      dpd_buf4_copy(&T2, CC_LAMPS, "LIJAB");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
      dpd_buf4_copy(&T2, CC_LAMPS, "Lijab");
      dpd_buf4_close(&T2);
  
      dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      dpd_buf4_copy(&T2, CC_LAMPS, "LIjAb");
      dpd_buf4_close(&T2);
    }
  }

  if (!params.ground) {
      /* multiply by R0 and create nonsymmetric L files */
      dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
      dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
      dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");

      dpd_file2_scm(&LIA, params.R0);
      dpd_file2_scm(&Lia, params.R0);
      dpd_buf4_scm(&LIJAB, params.R0);
      dpd_buf4_scm(&Lijab, params.R0);
      dpd_buf4_scm(&LIjAb, params.R0);

      /* add R1 and R2 */
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "RIA");
      dpd_file2_axpy(&R1, &LIA, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "Ria");
      dpd_file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "RIJAB");
      dpd_buf4_axpy(&R2, &LIJAB, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "Rijab");
      dpd_buf4_axpy(&R2, &Lijab, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, "RIjAb");
      dpd_buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_buf4_close(&R2);

      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "RIA");
      norm = dpd_file2_dot(&LIA, &R1);
      dpd_file2_close(&R1);

      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "Ria");
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);

      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "RIJAB");
      norm += dpd_buf4_dot(&LIJAB, &R2);
      dpd_buf4_close(&R2);

      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "Rijab");
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);

      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, "RIjAb");
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);

      fprintf(outfile,"initial norm <R|L> %15.10lf\n",norm);

      dpd_file2_scm(&LIA, 1.0/norm);
      dpd_file2_scm(&Lia, 1.0/norm);
      dpd_buf4_scm(&LIJAB, 1.0/norm);
      dpd_buf4_scm(&Lijab, 1.0/norm);
      dpd_buf4_scm(&LIjAb, 1.0/norm);

      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "RIA");
      norm = dpd_file2_dot(&LIA, &R1);
      dpd_file2_close(&R1);
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, "Ria");
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "RIJAB");
      norm += dpd_buf4_dot(&LIJAB, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, "Rijab");
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, "RIjAb");
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
      fprintf(outfile,"checking initial <R|L> normalization %15.10lf\n", norm);

      dpd_file2_close(&LIA);
      dpd_file2_close(&Lia);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_close(&Lijab);
      dpd_buf4_close(&LIjAb);
  }
}
