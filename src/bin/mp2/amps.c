#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double amps(void) 
{
  dpdbuf4 tIJAB;
  dpdbuf4 tijab;
  dpdbuf4 tIjAb;
  dpdbuf4 D;
  dpdbuf4 dIJAB;
  dpdbuf4 dijab;
  dpdbuf4 dIjAb;

  if(params.ref == 0) { /** RHF **/
  
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "tIjAb");
    dpd_buf4_close(&D);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&dIjAb, &tIjAb);
    dpd_buf4_close(&dIjAb);
    dpd_buf4_close(&tIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    
  }
  else if(params.ref == 2) { /** UHF **/
      
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_copy(&D, CC_TAMPS, "tIJAB");
    dpd_buf4_close(&D);
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_dirprd(&dIJAB, &tIJAB);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&dIJAB);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_copy(&D, CC_TAMPS, "tijab");
    dpd_buf4_close(&D);
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    dpd_buf4_dirprd(&dIJAB, &tIJAB);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&dIJAB);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "tIjAb");
    dpd_buf4_close(&D);
    dpd_buf4_init(&dIJAB, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_dirprd(&dIJAB, &tIJAB);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&dIJAB);
  }
  
}

