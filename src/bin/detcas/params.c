/*
** PARAMS.C
**
** Gets/prints user specified input
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <stdlib.h>
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>
#include "globals.h"

/*
** get_parameters(): Function gets the program running parameters such
**   as convergence.  These are stored in the Parameters data structure.
*/
void get_parameters(void)
{
  int i, errcod;
  char line1[133];
   
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":DETCAS");


  errcod = ip_string("DERTYPE", &(Params.dertype),0);
  if(errcod == IPE_KEY_NOT_FOUND) {
    Params.dertype = (char *) malloc(sizeof(char)*5);
    strcpy(Params.dertype, "NONE");
  }

  if (strcmp(Params.dertype, "NONE")==0) {
    Params.rms_grad_convergence = 4;
    Params.energy_convergence = 7;
  }
  else {
    Params.rms_grad_convergence = 7;
    Params.energy_convergence = 11;
  }

  /* Params.print_lvl is set in detcas.cc */
  Params.print_mos = 0;
  Params.filter_ints = 0;  /* assume we need all for MCSCF */
  Params.oei_file = PSIF_OEI;  /* contains frozen core operator */
  Params.oei_erase = 0;
  Params.tei_file = PSIF_MO_TEI;
  Params.tei_erase = 0;
  Params.opdm_file = PSIF_MO_OPDM;
  Params.opdm_erase = 0;
  Params.tpdm_file = PSIF_MO_TPDM;
  Params.tpdm_erase = 0;
  Params.lag_file = PSIF_MO_LAG;
  Params.lag_erase = 0;

  Params.fzc = 0;              /* no sense freezing orbs in their SCF form */
  Params.fci = 1;              /* probably usually want CASSCF */

  Params.scale_grad = 1;       /* scale the orbital gradient? */
  Params.diis_start = 3;       /* iteration to turn on DIIS */
  Params.diis_freq  = 1;       /* how often to do a DIIS extrapolation */
  Params.diis_min_vecs = 2;
  Params.diis_max_vecs = 8;
  Params.scale_step = 1.0;

  errcod = ip_data("PRINT","%d",&(Params.print_lvl),0);
  errcod = ip_boolean("PRINT_MOS",&(Params.print_mos),0);
  errcod = ip_data("CONVERGENCE","%d",
                   &(Params.rms_grad_convergence),0);
  errcod = ip_data("ENERGY_CONVERGENCE","%d",
                   &(Params.energy_convergence),0);
  errcod = ip_data("OEI_FILE","%d",&(Params.oei_file),0);
  errcod = ip_boolean("OEI_ERASE",&(Params.oei_erase),0);
  errcod = ip_data("TEI_FILE","%d",&(Params.tei_file),0);
  errcod = ip_data("LAG_FILE","%d",&(Params.lag_file),0);
  errcod = ip_boolean("TEI_ERASE",&(Params.tei_erase),0);
  errcod = ip_boolean("FREEZE_CORE",&(Params.fzc),0);
  errcod = ip_boolean("FCI",&(Params.fci),0);
  errcod = ip_data("OPDM_FILE","%d",&(Params.opdm_file),0);
  errcod = ip_data("TPDM_FILE","%d",&(Params.tpdm_file),0);
  errcod = ip_boolean("SCALE_GRAD",&(Params.scale_grad),0);
  errcod = ip_data("DIIS_START","%d",&(Params.diis_start),0);
  errcod = ip_data("DIIS_FREQ","%d",&(Params.diis_freq),0);
  errcod = ip_data("DIIS_MIN_VECS","%d",&(Params.diis_min_vecs),0);
  errcod = ip_data("DIIS_MAX_VECS","%d",&(Params.diis_max_vecs),0);
  errcod = ip_data("SCALE_STEP","%lf",&(Params.scale_step),0);
  errcod = ip_string("HESSIAN",&(Params.hessian),0);
  if (errcod == IPE_KEY_NOT_FOUND) {
    Params.hessian = (char *) malloc(sizeof(char)*12);
    strcpy(Params.hessian, "APPROX_DIAG");
  }
  if (strcmp(Params.hessian, "FULL")==0) {
    strcpy(Params.hessian, "APPROX_DIAG");
    fprintf(outfile, "(detcas): FULL Hessian not yet available\n");
  } 
  if ((strcmp(Params.hessian, "FULL")!=0) && (strcmp(Params.hessian, "DIAG")!=0) &&
      (strcmp(Params.hessian, "APPROX_DIAG")!=0)) {
    fprintf(outfile, "(detcas): Unrecognized Hessian option %s\n", Params.hessian);
    exit(0);
  }
}


/*
** print_parameters(): Function prints the program's running parameters
**   found in the Parameters structure.
*/
void print_parameters(void)
{
  fprintf(outfile, "\n") ;
  fprintf(outfile, "PARAMETERS: \n") ;
  fprintf(outfile, "   PRINT         =   %6d      PRINT_MOS     =   %6s\n", 
      Params.print_lvl, Params.print_mos ? "yes" : "no");
  fprintf(outfile, "   CONVERGENCE   =   %6d      E CONVERG     =   %6d\n",
      Params.rms_grad_convergence, Params.energy_convergence);
  fprintf(outfile, "   FCI           =   %6s      FREEZE CORE   =   %6s\n", 
      Params.fci ? "yes" : "no", Params.fzc ? "yes" : "no") ;
  fprintf(outfile, "   OEI FILE      =   %6d      OEI ERASE     =   %6s\n", 
      Params.oei_file, Params.oei_erase ? "yes" : "no");
  fprintf(outfile, "   TEI FILE      =   %6d      TEI ERASE     =   %6s\n", 
      Params.tei_file, Params.tei_erase ? "yes" : "no");
  fprintf(outfile, "   OPDM FILE     =   %6d      OPDM ERASE    =   %6s\n", 
      Params.lag_file, Params.opdm_erase ? "yes" : "no");
  fprintf(outfile, "   TPDM FILE     =   %6d      TPDM ERASE    =   %6s\n", 
      Params.tpdm_file, Params.tpdm_erase ? "yes" : "no");
  fprintf(outfile, "   LAG FILE      =   %6d      LAG ERASE     =   %6s\n", 
      Params.lag_file, Params.lag_erase ? "yes" : "no");
  fprintf(outfile, "   DIIS START    =   %6d      DIIS FREQ     =   %6d\n", 
      Params.diis_start, Params.diis_freq);
  fprintf(outfile, "   DIIS MIN VECS =   %6d      DIIS MAX VECS =   %6d\n", 
      Params.diis_min_vecs, Params.diis_max_vecs);
  fprintf(outfile, "   SCALE STEP    =   %6.2E    HESSIAN       =   %-12s\n",
      Params.scale_step, Params.hessian);
  fprintf(outfile, "\n") ;
  fflush(outfile) ;
}


