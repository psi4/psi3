#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Lnorm(void)
{
  dpdfile2 R1, L1, LIA, Lia, RIA, Ria;
  dpdbuf4 R2, L2, LIJAB, Lijab, LIjAb, RIJAB, Rijab, RIjAb;
  double norm;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    /*
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "New Lia");
    dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    */
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");

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
    fprintf(outfile,"initial norm <L|R>     %15.10lf\n",norm);

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
    fprintf(outfile,"checking normalization %15.10lf\n", norm);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }
}
