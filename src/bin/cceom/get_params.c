#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, iconv, ref;
  char *cachetype = NULL;
  char *junk;

  errcod = ip_string("REFERENCE", &(junk),0);
  if(!strcmp(junk, "RHF")) ref = 0;
  else if(!strcmp(junk, "ROHF")) ref = 1;
  else if(!strcmp(junk, "UHF")) ref = 2;
  else { 
    printf("Invalid value of input keyword REFERENCE: %s\n", junk);
    exit(2); 
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    printf("Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	   ref, params.ref);
    exit(2);
  }
  
  fndcor(&(params.memory),infile,outfile);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.cachetype = 0;

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
  if(params.ref == 2) /* No LRU cacheing yet for UHF references */
    params.cachetype = 0;

  params.cachetype = 0;

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tReference wfn   =    %4s\n", junk);
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tCache Level     =    %1d\n", 
          params.cachelev);
  fprintf(outfile, "\tCache Type      =    %4s\n", 
          params.cachetype ? "LOW" : "LRU");
  fprintf(outfile, "\n");

  free(junk);
}

