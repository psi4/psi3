#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void setup(void)
{
  int h;
  int *occpi, *virtpi;
  dpdfile2 FIJ, Fij, FAB, Fab;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;

  dpd_file2_init(&FIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&Fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&FAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&Fab, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&FIJ);
  dpd_file2_mat_init(&Fij);
  dpd_file2_mat_init(&FAB);
  dpd_file2_mat_init(&Fab);
  dpd_file2_mat_rd(&FIJ);
  dpd_file2_mat_rd(&Fij);
  dpd_file2_mat_rd(&FAB);
  dpd_file2_mat_rd(&Fab);


  dpd_file2_mat_close(&FIJ);
  dpd_file2_mat_close(&Fij);
  dpd_file2_mat_close(&FAB);
  dpd_file2_mat_close(&Fab);
  dpd_file2_close(&FIJ);
  dpd_file2_close(&Fij);
  dpd_file2_close(&FAB);
  dpd_file2_close(&Fab);
}
