#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, iconv;
  char *cachetype = NULL;
  char *junk;

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

  params.cachetype = 1;
  errcod = ip_string("CACHETYPE", &(cachetype),0);
  if(cachetype != NULL && strlen(cachetype)) {
      if(!strcmp(cachetype,"LOW")) params.cachetype = 1;
      else if(!strcmp(cachetype,"LRU")) params.cachetype = 0;
      else {
          fprintf(outfile, "Error in input: invalid CACHETYPE = %s\n",
                  cachetype);
          exit(1);
        }
      free(cachetype);
    }

  errcod = ip_string("REFERENCE", &(junk),0);
  if(!strcmp(junk, "RHF")) params.ref = 0;
  else if(!strcmp(junk, "ROHF")) params.ref = 1;
  else { 
     printf("Invalid value of input keyword REFERENCE: %s\n", junk);
     exit(2); 
    }
  
  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tMaxiter         =   %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence     = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart         =     %s\n", 
          params.restart ? "Yes" : "No");
  fprintf(outfile, "\tAO Basis        =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tCache Level     =    %1d\n", 
          params.cachelev);
  fprintf(outfile, "\tCache Type      =    %4s\n", 
          params.cachetype ? "LOW" : "LRU");
  fprintf(outfile, "\tReference wfn   =    %4s\n", junk);
  fprintf(outfile, "\n");
}

