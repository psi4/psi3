#define EXTERN
#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <file30_params.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

void start_io()
{
  infile = fopen("./input.dat", "r");
  if (overwrite_output)
      outfile = fopen("./output.dat", "w+");
  else
      outfile = fopen("./output.dat", "a+");
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_cwk_add(":INPUT");
  tstart(outfile);
  /*--- Initialize new IO system ---*/
  psio_init();

  return;
}

void stop_io()
{
  tstop(outfile);
  ip_done();
  psio_done();
  fclose(outfile);
  fclose(infile);
}


