#include <stdio.h>
#include <dpd.h>
#include <ccfiles.h>
#define EXTERN
#include "globals.h"

/* spinad_amps(): For RHF references, build the AA and BB amplitudes from 
** the AB amplitudes.
**
** T2(IJ,AB) = T2(ij,ab) = T2(Ij,Ab) - T2(Ij,Ba)
*/

void spinad_amps(void)
{
  dpdfile2 T1;
  dpdbuf4 T2AB1, T2AB2, T2;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&T1, CC_OEI, "New tia");
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2AB1, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_copy(&T2AB1, CC_TMP0, "tIjAb");
    dpd_buf4_sort(&T2AB1, CC_TMP0, pqsr, 0, 5, "tIjBa");
    dpd_buf4_close(&T2AB1);

    dpd_buf4_init(&T2AB1, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&T2AB2, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjBa");
    dpd_buf4_axpy(&T2AB2, &T2AB1, -1.0);
    dpd_buf4_close(&T2AB2);
    dpd_buf4_close(&T2AB1);

    dpd_buf4_init(&T2AB1, CC_TMP0, 0, 2, 7, 0, 5, 0, "tIjAb");
    /*  dpd_buf4_print(&T2AB1, outfile, 1); */
    dpd_buf4_copy(&T2AB1, CC_TAMPS, "New tIJAB");
    dpd_buf4_copy(&T2AB1, CC_TAMPS, "New tijab");
    dpd_buf4_close(&T2AB1);
  }
}
