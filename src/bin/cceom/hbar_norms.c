#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void hbar_norms() {
  double tval;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ;

  fprintf(outfile,"\n");

  dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
  dpd_file2_init(&Fae, CC_OEI, H_IRR, 1, 1, "Fae");
  tval = dpd_file2_dot_self(&FAE);
  tval += dpd_file2_dot_self(&Fae);
  dpd_file2_close(&Fae);
  dpd_file2_close(&FAE);
  fprintf(outfile,"Fae   dot Fae   total %15.10lf\n", tval);

  dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
  dpd_file2_init(&Fmi, CC_OEI, H_IRR, 0, 0, "Fmi");
  tval = dpd_file2_dot_self(&FMI);
  tval += dpd_file2_dot_self(&Fmi); 
  dpd_file2_close(&Fmi);
  dpd_file2_close(&FMI);
  fprintf(outfile,"Fmi   dot Fmi   total %15.10lf\n", tval);

  dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
  dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
  tval = dpd_file2_dot_self(&FME);
  tval += dpd_file2_dot_self(&Fme);
  dpd_file2_close(&Fme);
  dpd_file2_close(&FME);
  fprintf(outfile,"Fme   dot Fme   total %15.10lf\n", tval);


  dpd_buf4_init(&WMBIJ, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "WMBIJ");
  tval = 2 * dpd_buf4_dot_self(&WMBIJ);
  dpd_buf4_close(&WMBIJ);
  fprintf(outfile,"WMBIJ dot WMBIJ total %15.10lf\n", tval);

  dpd_buf4_init(&Wmbij, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "Wmbij");
  tval = 2 * dpd_buf4_dot_self(&Wmbij);
  dpd_buf4_close(&Wmbij);
  fprintf(outfile,"Wmbij dot Wmbij total %15.10lf\n", tval);

  dpd_buf4_init(&WMbIj, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
  tval = dpd_buf4_dot_self(&WMbIj);
  dpd_buf4_close(&WMbIj);
  fprintf(outfile,"WMbIj dot WMbIj total %15.10lf\n", tval);

  dpd_buf4_init(&WmBiJ, CC_HBAR, H_IRR, 11, 0, 11, 0, 0, "WmBiJ (Bm,Ji)");
  tval = dpd_buf4_dot_self(&WmBiJ);
  dpd_buf4_close(&WmBiJ);
  fprintf(outfile,"WmBiJ dot WmBiJ total %15.10lf\n", tval);


  return;
}
