#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

void get_eom_params()
{
  int errcod, i, j, sym, iconv,exist;

  /*  ip_cwk_clear();
      ip_cwk_add(progid); */

  eom_params.max_iter = 80;
  errcod = ip_data("MAX_ITER","%d",&(eom_params.max_iter),0);

/*
  eom_params.num_roots = 1;
  errcod = ip_data("NUM_ROOTS","%d",&(eom_params.num_roots),0);
*/

  eom_params.rpi = (int *) malloc(moinfo.nirreps * sizeof(int));
  if (ip_exist("RPI",0)) {
     ip_count("RPI", &i, 0);
     if (i != moinfo.nirreps) {
        fprintf(outfile,"Dim. of rpi vector must be %d\n", moinfo.nirreps) ;
        exit(0);
     }
     for (i=0;i<moinfo.nirreps;++i)
        errcod = ip_data("RPI","%d",&(eom_params.rpi[i]),1,i);
  }
  else { fprintf(outfile,"Must have rpi vector in input.\n"); exit(0); } 

  eom_params.prop_root = 1;
  for (i=0;i<moinfo.nirreps;++i)
     if (eom_params.rpi[i] > 0) { eom_params.prop_root = eom_params.rpi[i]; break; }
  errcod = ip_data("PROP_ROOT","%d",&(eom_params.prop_root),0);


  eom_params.excitation_range = 3;
  errcod = ip_data("EXCITATION_RANGE","%d",&(eom_params.excitation_range),0);

  eom_params.print_singles = 0;
  errcod = ip_data("PRINT_SINGLES","%d",&(eom_params.print_singles),0);

  eom_params.vectors_per_root_SS = 5;
  errcod = ip_data("VECTORS_PER_ROOT_SS","%d",&(eom_params.vectors_per_root_SS),0);

  eom_params.vectors_per_root = 5;
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

  eom_params.schmidt_add_residual_tol = 1E-2;
  errcod = ip_data("SCHMIDT_ADD_RESIDUAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.schmidt_add_residual_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.max_iter_SS = 100;


  fprintf(outfile, "\n\tCCEOM parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tRoots sought per irrep =");
  for (i=0;i<moinfo.nirreps;++i) fprintf(outfile, "%3d", eom_params.rpi[i]);
  fprintf(outfile,"\n");
  fprintf(outfile, "\tRoot for properties        = %5d\n", eom_params.prop_root);
  fprintf(outfile, "\tMax. number of iterations  = %5d\n", eom_params.max_iter);
  fprintf(outfile, "\tVectors stored per root    = %5d\n", eom_params.vectors_per_root);
  fprintf(outfile, "\tPrint HbarSS iterations?   = %5d\n", eom_params.print_singles);
  fprintf(outfile, "\tEigenvalue tolerance       = %5.1e\n", eom_params.eval_tol);
  fprintf(outfile, "\tResidual vector tolerance  = %5.1e\n", eom_params.residual_tol);
  fprintf(outfile, "\tComplex tolerance          = %5.1e\n", eom_params.complex_tol);
  fprintf(outfile, "\n\n");

  sym = 0;
  for (i=0;i<moinfo.nirreps;++i)
    for (j=0;j<moinfo.openpi[i];++j)
      sym = sym ^ i;

  moinfo.sym = sym;

/*
  fprintf(outfile,"Symmetry of ground state: %d : %s\n",
                         moinfo.sym,moinfo.labels[moinfo.sym]);
*/

}

