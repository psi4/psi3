#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod;

  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
}

