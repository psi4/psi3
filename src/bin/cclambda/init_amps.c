#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"


void init_amps(void)
{
  struct oe_dpdfile T1;
  struct dpdbuf T2;

  /* Restart from previous amplitudes if we can/should */
  /* Need to adjust this for new I/O
  if(params.restart && flen(CC_LIA) && flen(CC_Lia) && flen(CC_LIJAB)
     && flen(CC_Lijab) && flen(CC_LIjAb)) return;
     */

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_copy(&T1, CC_OEI, "LIA", 0, outfile);
  dpd_oe_file_close(&T1);

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_copy(&T1, CC_OEI, "Lia", 0, outfile);
  dpd_oe_file_close(&T1);

  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_copy(&T2, CC_LAMPS, "LIJAB", 0, outfile);
  dpd_buf_close(&T2);

  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_copy(&T2, CC_LAMPS, "Lijab", 0, outfile);
  dpd_buf_close(&T2);

  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_copy(&T2, CC_LAMPS, "LIjAb", 0, outfile);
  dpd_buf_close(&T2);
}
