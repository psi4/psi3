#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

void distribute(void);

void sort_tei(void)
{
  double tolerance;
  struct dpdfile A, B, C, D, E, F;
  int print;

  tolerance = params.tolerance;
  print = params.print_lvl;

  distribute();

  dpd_file_init(&A, CC_AINTS, 0, 0, "A <ij|kl>", 0, outfile);
  file_build(&A, 90, tolerance, 1, 1, 1, 0, print, outfile);
  dpd_file_close(&A);

  dpd_file_init(&B, CC_BINTS, 5, 5, "B <ab|cd>", 0, outfile);
  file_build(&B, 91, tolerance, 1, 1, 1, 0, print, outfile);
  dpd_file_close(&B);

  dpd_file_init(&C, CC_CINTS, 10, 10, "C <ia|jb>", 0, outfile);
  file_build(&C, 92, tolerance, 1, 1, 0, 0, print, outfile);
  dpd_file_close(&C);

  dpd_file_init(&D, CC_DINTS, 0, 5, "D <ij|ab>", 0, outfile);
  file_build(&D, 93, tolerance, 0, 0, 1, 0, print, outfile);
  dpd_file_close(&D);

  dpd_file_init(&E, CC_EINTS, 11, 0, "E <ai|jk>", 0, outfile);
  file_build(&E, 94, tolerance, 0, 1, 0, 0, print, outfile);
  dpd_file_close(&E);

  dpd_file_init(&F, CC_FINTS, 10, 5, "F <ia|bc>", 0, outfile);
  file_build(&F, 95, tolerance, 0, 1, 0, 0, print, outfile);
  dpd_file_close(&F);

}
