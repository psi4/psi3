#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, iconv;

  params.maxiter = 50;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-12;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);
  params.restart = 1;
  errcod = ip_boolean("RESTART", &(params.restart),0);
  /* If the MO orbital phases are screwed up, don't restart */
  if(!moinfo.phase) params.restart = 0;

  fndcor(&(params.memory),infile,outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);
  params.aobasis = 0;  /* AO basis code not yet working for lambda */
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter     =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart     =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\tCache Level =     %1d\n", params.cachelev);
  fprintf(outfile, "\tAO Basis        =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\n");
}

