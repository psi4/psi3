#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Lsave(void)
{
  dpdfile2 L1;
  dpdbuf4 L2;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&L1, CC_OEI, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "New Lia");
    dpd_file2_copy(&L1, CC_OEI, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIjAb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&L1, CC_OEI, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_OEI, L_irr, 2, 3, "New Lia");
    dpd_file2_copy(&L1, CC_OEI, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIjAb");
    dpd_buf4_close(&L2);

  }
}

