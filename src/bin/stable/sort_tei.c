#include <stdio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void distribute(void);

int file_build(dpdfile4 *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs, int keep);

void sort_tei(void)
{
  double tolerance;
  dpdfile4 A, C, D, E;

  tolerance = params.tolerance;

  distribute();

  dpd_file4_init_nocache(&A, PSIF_MO_HESS, 0, 0, 0, "A <ij|kl>");
  file_build(&A, 90, tolerance, 1, 1, 1, 0);
  dpd_file4_close(&A);

  dpd_file4_init_nocache(&C, PSIF_MO_HESS, 0, 10, 10, "C <ia|jb>");
  file_build(&C, 91, tolerance, 1, 1, 0, 0);
  dpd_file4_close(&C);

  dpd_file4_init_nocache(&D, PSIF_MO_HESS, 0, 0, 5, "D <ij|ab>");
  file_build(&D, 92, tolerance, 0, 0, 1, 0);
  dpd_file4_close(&D);

  dpd_file4_init_nocache(&E, PSIF_MO_HESS, 0, 11, 0, "E <ai|jk>");
  file_build(&E, 93, tolerance, 0, 1, 0, 0);
  dpd_file4_close(&E);

  fflush(outfile);

}
