#define EXTERN
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <stdlib.h>
#include <file30_params.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

void start_io()
{
  infile = fopen("./input.dat", "r");
  outfile = fopen("./output.dat", "w");
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":INPUT");
  ip_cwk_add(":DEFAULT");
  tstart(outfile);
  /*--- Initialize new IO system ---*/
  /* psio_init(); */

  return;
}

void stop_io()
{
  tstop(outfile);
  ip_done();
  /* psio_done();*/
  fclose(outfile);
  fclose(infile);
}


