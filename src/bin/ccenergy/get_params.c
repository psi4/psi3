#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, iconv, forceit;
  char *cachetype = NULL;
  char *junk;
  int *mu_irreps;

  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "CCSD") && strcmp(params.wfn, "CCSD_T") &&
     strcmp(params.wfn, "EOM_CCSD") && strcmp(params.wfn, "LEOM_CCSD") &&
     strcmp(params.wfn, "BCCD") && strcmp(params.wfn,"BCCD_T") &&
     strcmp(params.wfn,"CC3")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  if(!strcmp(params.wfn,"BCCD") || !strcmp(params.wfn,"BCCD_T")) 
    params.brueckner = 1;
  else params.brueckner = 0;

  params.semicanonical = 0;
  errcod = ip_string("REFERENCE", &(junk),0);
  /* if no reference is given, assume rhf */
  if (errcod != IPE_OK) {
    params.ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) params.ref = 0;
    else if(!strcmp(junk,"ROHF") && !strcmp(params.wfn,"MP2") || !strcmp(params.wfn,"CCSD_T")) {
      params.ref = 2;
      params.semicanonical = 1;
    }
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

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  params.maxiter = 50;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-7;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);
  params.restart = 1;
  errcod = ip_boolean("RESTART", &(params.restart),0);
  /* If the MO orbital phases are screwed up, don't restart */
  if(!moinfo.phase) params.restart = 0;
  /* BUT, the user can force an override of the phase problem */
  forceit = 0;
  errcod = ip_boolean("FORCE_RESTART", &forceit,0);
  if(forceit) params.restart = 1;

  fndcor(&(params.memory),infile,outfile);

  if(ip_exist("AO_BASIS",0)) {
    errcod = ip_string("AO_BASIS", &(params.aobasis),0);
  }
  else params.aobasis = strdup("NONE");
  if(strcmp(params.aobasis,"DISK") && strcmp(params.aobasis,"DIRECT") &&
     strcmp(params.aobasis,"NONE")) {
    fprintf(outfile, "Error in input: invalid AO_BASIS = %s\n",
	    params.aobasis);
    exit(PSI_RETURN_FAILURE);
  }

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
      exit(PSI_RETURN_FAILURE);
    }
    free(cachetype);
  }
  if(params.ref == 2) /* No LOW cacheing yet for UHF references */
    params.cachetype = 0;

  params.diis = 1;
  errcod = ip_boolean("DIIS", &(params.diis),0);

  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);
  local.cutoff = 0.02;
  errcod = ip_data("LOCAL_CUTOFF", "%lf", &(local.cutoff), 0);

  if(ip_exist("LOCAL_METHOD",0)) {
    errcod = ip_string("LOCAL_METHOD", &(local.method), 0);
    if(strcmp(local.method,"AOBASIS") && strcmp(local.method,"WERNER")) {
      fprintf(outfile, "Invalid local correlation method: %s\n", local.method);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local) {
    local.method = (char *) malloc(7 * sizeof(char));
    sprintf(local.method, "%s", "WERNER");
  }

  if(ip_exist("LOCAL_WEAKP",0)) {
    errcod = ip_string("LOCAL_WEAKP", &(local.weakp), 0);
    if(strcmp(local.weakp,"MP2") && strcmp(local.weakp,"NEGLECT") && strcmp(local.weakp,"NONE")) {
      fprintf(outfile, "Invalid method for treating local pairs: %s\n", local.weakp);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local) {
    local.weakp = (char *) malloc(4 * sizeof(char));
    sprintf(local.weakp, "%s", "MP2");
  }

  local.filter_singles = 1;
  ip_boolean("LOCAL_FILTER_SINGLES", &(local.filter_singles), 0);

  local.cphf_cutoff = 0.01;
  ip_data("LOCAL_CPHF_CUTOFF", "%lf", &(local.cphf_cutoff), 0);

  local.freeze_core = NULL;
  ip_string("FREEZE_CORE", &local.freeze_core, 0);
  if(local.freeze_core == NULL) local.freeze_core = strdup("FALSE");

  if(ip_exist("LOCAL_PAIRDEF",0)){
    errcod = ip_string("LOCAL_PAIRDEF", &(local.pairdef), 0);
    if(strcmp(local.pairdef,"BP") && strcmp(local.pairdef,"RESPONSE")) {
      fprintf(outfile, "Invalid keyword for strong/weak pair definition: %s\n", local.pairdef);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local)
    local.pairdef = strdup("RESPONSE");

  params.num_amps = 10;
  if(ip_exist("NUM_AMPS",0)) {
    errcod = ip_data("NUM_AMPS", "%d", &(params.num_amps), 0);
  }

  params.bconv = 1e-5;
  errcod = ip_data("BRUECKNER_CONV", "%d", &(iconv), 0);
  if(errcod == IPE_OK) params.bconv = 1.0*pow(10.0,(double) -iconv);

  params.print_mp2_amps = 0;
  errcod = ip_boolean("PRINT_MP2_AMPS", &(params.print_mp2_amps), 0);

  params.analyze = 0;
  errcod = ip_boolean("ANALYZE", &(params.analyze), 0);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function   =    %6s\n", params.wfn);
  if(params.semicanonical) {
    fprintf(outfile, "\tReference wfn   =    ROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
    fprintf(outfile, "\tReference wfn   =    %5s\n",
	    (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  }
  if(params.brueckner) 
    fprintf(outfile, "\tBrueckner conv. =    %3.1e\n", params.bconv);
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tMaxiter         =   %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence     = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart         =     %s\n", 
	  params.restart ? "Yes" : "No");
  fprintf(outfile, "\tDIIS            =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff       = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method      =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs        =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef);
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }
  fprintf(outfile, "\tAO Basis        =     %s\n", params.aobasis);
  fprintf(outfile, "\tCache Level     =    %1d\n", params.cachelev);
  fprintf(outfile, "\tCache Type      =    %4s\n", 
	  params.cachetype ? "LOW" : "LRU");
  fprintf(outfile, "\tPrint Level     =    %1d\n",  params.print);
  fprintf(outfile, "\t# Amps to Print =    %1d\n",  params.num_amps);
  fprintf(outfile, "\tPrint MP2 Amps? =    %s\n",  params.print_mp2_amps ?
	  "Yes" : "No" );
  fprintf(outfile, "\tAnalyze T2 Amps =    %s\n",  params.analyze ? "Yes" : "No" );
  fprintf(outfile, "\n");

}

