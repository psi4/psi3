#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"


void init_amps(void)
{
  dpdfile2 T1;
  dpdbuf4 T2;

  /* Restart from previous amplitudes if we can/should */
  /* Need to adjust this for new I/O
     if(params.restart && flen(CC_LIA) && flen(CC_Lia) && flen(CC_LIJAB)
     && flen(CC_Lijab) && flen(CC_LIjAb)) return;
  */

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
