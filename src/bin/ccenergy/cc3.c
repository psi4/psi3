#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

void cc3_RHF(void);
void cc3_UHF_AAA(void);
void cc3_UHF_BBB(void);
void cc3_UHF_AAB(void);
void cc3_UHF_BBA(void);

void cc3(void)
{
  dpdfile2 t1;
  dpdbuf4 T2;

  if(params.ref == 0) cc3_RHF();
  else if(params.ref == 2) {
    cc3_UHF_AAA();
    cc3_UHF_BBB();
    cc3_UHF_AAB();
    cc3_UHF_BBA();

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2);

  }
}
