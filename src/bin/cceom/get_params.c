#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, iconv;
  char *cachetype = NULL;
  char *read_ref, *read_eom_ref;

  errcod = ip_string("REFERENCE", &(read_ref),0);
  if(!strcmp(read_ref, "RHF")) params.ref = 0;
  else if(!strcmp(read_ref, "ROHF")) params.ref = 1;
  else if(!strcmp(read_ref, "UHF")) params.ref = 2;
  else { 
    fprintf(outfile,
        "\nInvalid value of input keyword REFERENCE: %s\n", read_ref);
    exit(2); 
  }

  if (params.ref == 0) { // for RHF refs, allow CCEOM to do RHF, ROHF, UHF modes
    errcod = ip_string("EOM_REFERENCE", &(read_eom_ref),0);
    if (errcod == IPE_OK) {
      if(!strcmp(read_eom_ref, "RHF")) params.eom_ref = 0;
      else if(!strcmp(read_eom_ref, "ROHF")) params.eom_ref = 1;
      else if(!strcmp(read_eom_ref, "UHF")) params.eom_ref = 2;
      else { 
        fprintf(outfile,
            "\nInvalid value of input keyword EOM_REFERENCE: %s\n", read_eom_ref);
        exit(2); 
      }
    }
    else {
      params.eom_ref = 0;
      read_eom_ref = (char *) malloc(10*sizeof(char));
      sprintf(read_eom_ref,"%s","RHF"); // just for printing below
    }
  }
  else if (params.ref == 1) { // for ROHF refs, allow CCEOM to do ROHF & UHF modes
    errcod = ip_string("EOM_REFERENCE", &(read_eom_ref),0);
    if (errcod == IPE_OK) {
      if(!strcmp(read_eom_ref, "ROHF")) params.eom_ref = 1;
      else if(!strcmp(read_eom_ref, "UHF")) params.eom_ref = 2;
      else { 
        fprintf(outfile,
            "\nInvalid value of input keyword EOM_REFERENCE: %s\n", read_eom_ref);
        exit(2); 
      }
    }
    else {
      params.eom_ref = 1;
      read_eom_ref = (char *) malloc(10*sizeof(char));
      sprintf(read_eom_ref,"%s","ROHF"); // just for printing below
    }
  }
  else { // run in UHF mode - ignore EOM_REFERENCE
    params.eom_ref = 2;
    read_eom_ref = (char *) malloc(10*sizeof(char));
    sprintf(read_eom_ref,"%s","UHF"); // just for printing below
  }

  if (params.eom_ref == 2) {
    fprintf(outfile, "\nCCEOM not yet UHF capable\n");
    exit(2); 
  }

  /* if(params.ref != ref) {
     fprintf(outfile,
     "\nValue of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
     ref, params.ref);
     } */

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

  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);

  local.cutoff = 1e-2;
  errcod = ip_data("LOCAL_CUTOFF", "%d", &(iconv), 0);
  if(errcod == IPE_OK) local.cutoff = 1.0 * pow(10.0, (double) -iconv);

  if(ip_exist("LOCAL_METHOD",0)) {
    errcod = ip_string("LOCAL_METHOD", &(local.method), 0);
    if(strcmp(local.method,"AOBASIS") && strcmp(local.method,"WERNER")) {
      fprintf(outfile, "\nInvalid local correlation method: %s\n", local.method);
      exit(2);
    }
  }
  else if(params.local) {
    local.method = (char *) malloc(7 * sizeof(char));
    sprintf(local.method, "%s", "WERNER");
  }

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tReference wfn   =    %4s\n", read_ref);
  fprintf(outfile, "\tReference EOM wfn=    %4s\n", read_eom_ref);
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =     %s\n", 
      params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tCache Level     =    %1d\n", 
      params.cachelev);
  fprintf(outfile, "\tCache Type      =    %4s\n", 
      params.cachetype ? "LOW" : "LRU");
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
  }
  fprintf(outfile, "\n");

  free(read_ref);
  free(read_eom_ref);
}

