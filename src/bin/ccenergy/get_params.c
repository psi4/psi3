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
  params.convergence = 1e-8;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);
  params.restart = 1;
  errcod = ip_boolean("RESTART", &(params.restart),0);
  /* If the MO orbital phases are screwed up, don't restart */
  if(!moinfo.phase) params.restart = 0;

  fndcor(&(params.memory),infile,outfile);

  params.aobasis = 0;
  errcod = ip_boolean("AOBASIS", &(params.aobasis),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter     =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart     =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\tAO Basis    =     %s\n", params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tCache Level =    %1d\n", params.cachelev);
  fprintf(outfile, "\n");
}

