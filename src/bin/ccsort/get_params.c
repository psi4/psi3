#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void get_params()
{
  int i, errcod, tol;
  char *junk;

  errcod = ip_string("REFERENCE", &(junk),0);
  if(!strcmp(junk, "RHF")) params.ref = 0;
  else if(!strcmp(junk, "ROHF")) params.ref = 1;
  else if(!strcmp(junk, "UHF")) params.ref = 2;
  else { 
    printf("Invalid value of input keyword REFERENCE: %s\n", junk);
    exit(2); 
  }
  free(junk);

  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) params.dertype = 0;
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(2); 
    }
    free(junk);
  }

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);

  params.print_lvl = 1;
  errcod = ip_data("PRINT_LVL","%d",&(params.print_lvl),0);

  params.keep_TEIFile = 1;
  errcod = ip_boolean("KEEP_TEIFILE",&(params.keep_TEIFile),0);

  params.keep_OEIFile = 1;
  errcod = ip_boolean("KEEP_OEIFILE",&(params.keep_OEIFile),0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE", "%d", &(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tReference wfn   =    %5s\n", 
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tDerivative      =    %5s\n", 
	  (params.dertype == 0) ? "None" : "First");
  fprintf(outfile, "\tMemory (Mbytes) =    %5.1f\n", params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =    %5s\n", params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tCache Level     =    %5d\n", params.cachelev);
  fprintf(outfile, "\tCache Type      =    %5s\n", "LRU");
  fprintf(outfile, "\n");
}
