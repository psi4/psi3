#include <stdio.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
}

