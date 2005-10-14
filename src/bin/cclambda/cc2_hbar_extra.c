#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void cc2_hbar_extra(void) {
  dpdbuf4 W1, W2, W;

  if(!strcmp(params.wfn,"CC2")) {
    if(params.ref == 0) { /** RHF **/
      /* 2 W(ME,jb) + W(Me,Jb) */
      dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
      dpd_buf4_copy(&W1, CC2_HET1, "CC2 2 W(ME,jb) + W(Me,Jb)");
      dpd_buf4_close(&W1);
      dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
      dpd_buf4_init(&W2, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
      dpd_buf4_axpy(&W2, &W1, 2);
      dpd_buf4_close(&W2);
      dpd_buf4_close(&W1);
    }
  }
}
