#include <libdpd/dpd.h>
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#define EXTERN
#include "globals.h"

void energy(void)
{
  double energy;
  dpdbuf4 I, D, T2A;
  
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&I, CC_TMP0, "tIjAb");
  dpd_buf4_close(&I);

  dpd_buf4_init(&T2A, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D, &T2A);
  dpd_buf4_close(&D);

  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  energy = dpd_buf4_dot(&I, &T2A);
  dpd_buf4_close(&T2A);
  dpd_buf4_close(&I);

  fprintf(outfile,"\n");
  fprintf(outfile,"\tMP2 correlation energy \t=\t %.12f\n",energy);
  fprintf(outfile,"\tMP2 total energy       \t=\t%.12f\n",mo.Escf+energy);
}
