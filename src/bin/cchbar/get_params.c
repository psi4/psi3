#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod;

  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  errcod = ip_string("WFN", &(params.wfn), 0);

}

