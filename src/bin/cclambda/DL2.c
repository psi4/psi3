#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void DL2(void)
{
  struct dpdbuf D;

  /* RHS = <ij||ab> */
  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
               0, outfile);
  dpd_copy(&D, CC_LAMPS, "New LIJAB",0,outfile);
  dpd_copy(&D, CC_LAMPS, "New Lijab",0,outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&D, CC_DINTS,0,5,0,5,0,"D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_LAMPS, "New LIjAb",0,outfile);
  dpd_buf_close(&D);

}
