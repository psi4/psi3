#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Lsave(void)
{
  struct oe_dpdfile L1;
  struct dpdbuf L2;

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "New LIA", 0, outfile);
  dpd_oe_copy(&L1, CC_OEI, "LIA", 0, outfile);
  dpd_oe_file_close(&L1);

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "New Lia", 0, outfile);
  dpd_oe_copy(&L1, CC_OEI, "Lia", 0, outfile);
  dpd_oe_file_close(&L1);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "New LIJAB", 0, outfile);
  dpd_copy(&L2, CC_LAMPS, "LIJAB", 0, outfile);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "New Lijab", 0, outfile);
  dpd_copy(&L2, CC_LAMPS, "Lijab", 0, outfile);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb", 0, outfile);
  dpd_copy(&L2, CC_LAMPS, "LIjAb", 0, outfile);
  dpd_buf_close(&L2);
}

