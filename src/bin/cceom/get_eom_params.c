#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_eom_params()
{
  int errcod, i, j, sym, iconv,exist, state_irrep, c_irrep;

  eom_params.max_iter = 80;
  errcod = ip_data("MAX_ITER","%d",&(eom_params.max_iter),0);

  eom_params.states_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
  if (ip_exist("STATES_PER_IRREP",0)) {
     ip_count("STATES_PER_IRREP", &i, 0);
     if (i != moinfo.nirreps) {
        fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps) ;
        exit(0);
     }
     for (i=0;i<moinfo.nirreps;++i)
        errcod = ip_data("STATES_PER_IRREP","%d",&(eom_params.states_per_irrep[i]),1,i);
  }
  else { fprintf(outfile,"Must have states_per_irrep vector in input.\n"); exit(0); } 

  eom_params.cs_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
  for (state_irrep=0; state_irrep<moinfo.nirreps; ++state_irrep) {
    for (c_irrep=0;c_irrep<moinfo.nirreps;++c_irrep)
       if ((moinfo.sym ^ c_irrep) == state_irrep) 
           eom_params.cs_per_irrep[c_irrep] = eom_params.states_per_irrep[state_irrep];
  }

  /* eom_params.prop_sym holds irrep of state used for properties */
  if (ip_exist("PROP_SYM",0)) {
    ip_data("PROP_SYM","%d",&(eom_params.prop_sym),0);
    eom_params.prop_sym = eom_params.prop_sym - 1;
  }
  else { /* assume last irrep is right one */
    for (i=0;i<moinfo.nirreps;++i)
      if (eom_params.states_per_irrep[i]) eom_params.prop_sym = i;
  }

  if (ip_exist("PROP_ROOT",0)) {
    ip_data("PROP_ROOT","%d",&(eom_params.prop_root),0);
    if (eom_params.prop_root > eom_params.states_per_irrep[eom_params.prop_sym]) {
      fprintf(outfile,"prop_root is too big\n");
      exit(1);
    }
  }
  else {
    eom_params.prop_root = eom_params.states_per_irrep[eom_params.prop_sym];
  }

  eom_params.excitation_range = 2;
  errcod = ip_data("EXCITATION_RANGE","%d",&(eom_params.excitation_range),0);

  eom_params.print_singles = 0;
  errcod = ip_data("PRINT_SINGLES","%d",&(eom_params.print_singles),0);

  eom_params.vectors_per_root_SS = 5;
  errcod = ip_data("VECTORS_PER_ROOT_SS","%d",&(eom_params.vectors_per_root_SS),0);

  eom_params.vectors_per_root = 6;
  errcod = ip_data("VECTORS_PER_ROOT","%d",&(eom_params.vectors_per_root),0);

  eom_params.complex_tol = 1E-12;
  errcod = ip_data("COMPLEX_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.complex_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.residual_tol = 1E-4;
  errcod = ip_data("RESIDUAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.residual_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.residual_tol_SS = 1E-6;
  errcod = ip_data("RESIDUAL_TOL_SS","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.residual_tol_SS = 1.0*pow(10.0,(double) -iconv);

  eom_params.eval_tol = 1E-6;
  errcod = ip_data("EVAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.eval_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.schmidt_add_residual_tol = 1E-3;
  errcod = ip_data("SCHMIDT_ADD_RESIDUAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.schmidt_add_residual_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.max_iter_SS = 100;


  fprintf(outfile, "\n\tCCEOM parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tStates sought per irrep     =");
  for (i=0;i<moinfo.nirreps;++i) fprintf(outfile, " %s %d,", moinfo.labels[i],
      eom_params.states_per_irrep[i]);
  fprintf(outfile,"\n");
  fprintf(outfile, "\tMax. number of iterations   = %5d\n", eom_params.max_iter);
  fprintf(outfile, "\tVectors stored per root     = %5d\n", eom_params.vectors_per_root);
  fprintf(outfile, "\tPrint HbarSS iterations?    = %5d\n", eom_params.print_singles);
  fprintf(outfile, "\tExcitation range for HBarSS = %5d\n", eom_params.excitation_range);
  fprintf(outfile, "\tEigenvalue tolerance        = %5.1e\n", eom_params.eval_tol);
  fprintf(outfile, "\tResidual vector tolerance   = %5.1e\n", eom_params.residual_tol);
  fprintf(outfile, "\tComplex tolerance           = %5.1e\n", eom_params.complex_tol);
  fprintf(outfile, "\tRoot for properties         = %5d\n", eom_params.prop_root);
  fprintf(outfile, "\tSym of state for properties = %6s\n", moinfo.labels[eom_params.prop_sym]);
  fprintf(outfile, "\n\n");

}

