#include <stdio.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, tol;

  params.tpdmfile = PSIF_MO_TPDM;
  errcod = ip_data("TPDM_FILE", "%d", &(params.tpdmfile),0);
  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE","%d",&(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory),infile,outfile);
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tTwo-pdm file =    %4d\n", params.tpdmfile);
  fprintf(outfile, "\tTolerance    = %3.1e\n", params.tolerance);
  fprintf(outfile, "\n");
}

