#include <stdio.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
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
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter     =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart     =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\n");
}

