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
  char *junk;

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

  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) params.dertype = 0;
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }

  if ( (!strcmp(params.wfn,"EOM_CCSD")) && (params.dertype == 0) )
    params.connect_xi = 0;
  else
    params.connect_xi = 1;
  errcod = ip_boolean("CONNECT_XI",&(params.connect_xi),0);
  
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
  fprintf(outfile, "\tXi connected=     %s\n", 
          (params.connect_xi) ? "Yes" : "No");
  fprintf(outfile, "\n");
}

