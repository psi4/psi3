#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_params(void)
{
  int errcod, iconv,i;
  char lbl[32];

  params.maxiter = 50;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-7;
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

  /* determine number of left-hand states per irrep */
  params.states_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
  for (i=0;i<moinfo.nirreps;++i) params.states_per_irrep[i] = 0;

  if (params.ground) { /* for ground state only find 1 A1 lambda */
    params.states_per_irrep[moinfo.sym] = 1;
  }
  else { /* excited state */
    if (ip_exist("STATES_PER_IRREP",0)) {
      ip_count("STATES_PER_IRREP", &i, 0);
      if (i != moinfo.nirreps) {
        fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps) ;
        exit(0);
      }
      for (i=0;i<moinfo.nirreps;++i)
        errcod = ip_data("STATES_PER_IRREP","%d",&(params.states_per_irrep[i]),1,i);
    }
    else { fprintf(outfile,"Must have states_per_irrep vector in input.\n"); exit(0); }
  }

  /* determine Ls per irrep */
  params.Ls_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
  for (i=0;i<moinfo.nirreps;++i) {
    params.Ls_per_irrep[i^moinfo.sym] =  params.states_per_irrep[i];
  }

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

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter     =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart     =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\tCache Level =     %1d\n", params.cachelev);
  fprintf(outfile, "\tAO Basis        =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tExcited State Computation =     %s\n", 
          params.ground ? "No" : "Yes");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
  }
  fprintf(outfile, "\tStates sought per irrep     =");
  for (i=0;i<moinfo.nirreps;++i)
    fprintf(outfile, " %s %d,", moinfo.labels[i],
        params.states_per_irrep[i]);
  fprintf(outfile,"\n");
  fprintf(outfile, "\tLs sought per irrep         =");
  for (i=0;i<moinfo.nirreps;++i)
    fprintf(outfile, " %s %d,", moinfo.labels[i],
        params.Ls_per_irrep[i]);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  fprintf(outfile, "\n");

  /* compute total number of states */
  params.nstates = 0;
  for (i=0;i<moinfo.nirreps;++i)
    params.nstates += params.states_per_irrep[i]; 

  /* get cceom energies for all states from CC_INFO */
  params.cceom_energy = init_array(params.nstates+1);
  params.R0 = init_array(params.nstates+1);
  if (!params.ground) {
    params.cceom_energy[0] = 0.0; /* ground state */
    params.R0[0] = 1.0;
    for (i=1; i<=params.nstates; ++i) {
      sprintf(lbl,"EOM CCSD Energy for root %d", i);
      psio_read_entry(CC_INFO, lbl, (char *) &(params.cceom_energy[i]),sizeof(double));
      sprintf(lbl,"EOM CCSD R0 for root %d", i);
      psio_read_entry(CC_INFO, lbl, (char *) &(params.R0[i]),sizeof(double));
    }
  }
  else { /* just a ground state calculation */
    params.cceom_energy[0] = 0.0;
    params.R0[0] = 1;
  }

#ifdef EOM_DEBUG
  if (!params.ground) {
    fprintf(outfile,"\tNstates = %d\n", params.nstates);
    for (i=1; i<=params.nstates; ++i) {
      fprintf(outfile,"\tEnergy and R0 for root %d: %15.10lf %15.10lf\n",
        i, params.cceom_energy[i], params.R0[i]);
    }
  }
#endif

  /* determine L0 value */
  params.L0 = (params.ground ? 1.0 : 0.0);

  return;
}

