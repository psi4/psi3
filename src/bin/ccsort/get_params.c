#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  int errcod, tol;
  char *junk;
  
  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "MP2") && strcmp(params.wfn, "CCSD") && 
     strcmp(params.wfn, "CCSD_T") && strcmp(params.wfn, "EOM_CCSD") && 
     strcmp(params.wfn, "LEOM_CCSD") && strcmp(params.wfn, "BCCD") && 
     strcmp(params.wfn,"BCCD_T") && strcmp(params.wfn, "SCF") &&
     strcmp(params.wfn,"CIS") && strcmp(params.wfn,"RPA")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  /* NB: SCF wfns are allowed because, at present, ccsort is needed for 
     RPA calculations */
  
  errcod = ip_string("REFERENCE", &(junk),0);
  if (errcod != IPE_OK) {
    /* if no reference is given, assume rhf */
    params.ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) params.ref = 0;
    else if(!strcmp(junk, "ROHF")) params.ref = 1;
    else if(!strcmp(junk, "UHF")) params.ref = 2;
    else { 
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }

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

  if(ip_exist("AO_BASIS",0)) {
    errcod = ip_string("AO_BASIS", &(params.aobasis),0);
  }
  else params.aobasis = strdup("NONE");

  if(!strcmp(params.wfn,"MP2") || !strcmp(params.aobasis,"DISK") ||
     !strcmp(params.aobasis,"DIRECT")) {
    params.make_abcd = 0;
  }
  else {
    params.make_abcd = 1;
  }
  
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
  fprintf(outfile, "\tWave function   =\t%s\n", params.wfn);
  fprintf(outfile, "\tReference wfn   =\t%s\n", 
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  if(params.dertype == 0) fprintf(outfile, "\tDerivative      =\tNone\n");
  else if(params.dertype == 1) fprintf(outfile, "\tDerivative      =\tFirst\n");
  else if(params.dertype == 3) fprintf(outfile, "\tDerivative      =\tResponse\n");
  fprintf(outfile, "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =\t%s\n", params.aobasis);
  fprintf(outfile, "\tMake (ab|cd)    =\t%s\n", 
	  (params.make_abcd == 1) ? "True" : "False");
  fprintf(outfile, "\tCache Level     =\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type      =\t%s\n", "LRU");
  fprintf(outfile, "\n");
}
