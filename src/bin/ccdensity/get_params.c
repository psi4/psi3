#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, tol;

  errcod = ip_string("WFN", &(params.wfn), 0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE","%d",&(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);
  
  params.relax_opdm = 1;
  errcod = ip_boolean("RELAX_OPDM", &(params.relax_opdm),0);
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tTolerance    = %3.1e\n", params.tolerance);
  fprintf(outfile, "\tCache Level  =    %1d\n", params.cachelev);
  fprintf(outfile, "\tAO Basis     =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tRelax OPDM   =     %s\n", 
          params.relax_opdm ? "Yes" : "No");
  fprintf(outfile, "\tExcited State=     %s\n", 
          (!params.ground) ? "Yes" : "No");
  fprintf(outfile, "\tCompute Xi   =     %s\n", 
          (params.calc_xi) ? "Yes" : "No");
  fprintf(outfile, "\tUse Zeta    =     %s\n", 
          (params.use_zeta) ? "Yes" : "No");
  fprintf(outfile, "\n");
}

