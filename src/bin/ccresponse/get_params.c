#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod, ref, count, iconv, *tmpi;
  char *junk, units[20];

  params.print = 1;
  errcod = ip_data("PRINT","%d",&(params.print),0);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  params.cachelev = 0;

  errcod = ip_string("REFERENCE", &(junk),0);
  /* if no reference is given, assume rhf */
  if (errcod != IPE_OK) {
    ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) ref = 0;
    else if(!strcmp(junk, "ROHF")) ref = 1;
    else if(!strcmp(junk, "UHF")) ref = 2;
    else { 
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	    ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
    params.ref = ref;
  }

  /* grab the field strength from input -- a few different units are converted to E_h */
  params.omega = 0.0; /* static polarizability by default */
  if(ip_exist("OMEGA",0)) {
    errcod = ip_count("OMEGA", &count, 0);

    if(errcod == IPE_NOT_AN_ARRAY)  /* assume Hartrees */
      errcod = ip_data("OMEGA", "%lf", &(params.omega), 0);

    else if(count == 2) {

      errcod = ip_data("OMEGA", "%lf", &(params.omega), 1, 0);
      errcod = ip_data("OMEGA", "%s", units, 1, 1);

      for(junk = units; *junk != '\0'; junk++)
	if(*junk>='a' && *junk <= 'z') *junk += 'A' - 'a';

      if(!strcmp(units, "HZ")) params.omega *= _h / _hartree2J;
      else if(!strcmp(units, "NM")) params.omega = (_c*_h*1e9)/(params.omega*_hartree2J);
      else if(!strcmp(units, "EV")) params.omega /= _hartree2ev;
    }
  }

  moinfo.mu_irreps = init_int_array(3);
  errcod = ip_int_array("MU_IRREPS", moinfo.mu_irreps, 3);
  if(errcod == IPE_OK) {
    moinfo.irrep_x = moinfo.mu_irreps[0];
    moinfo.irrep_y = moinfo.mu_irreps[1];
    moinfo.irrep_z = moinfo.mu_irreps[2];
  }
  else {
    fprintf(outfile, "\nYou must supply the irreps of x, y, and z with the MU_IRREPS keyword.\n");
    exit(PSI_RETURN_FAILURE);
  }

  moinfo.l_irreps = init_int_array(3);
  errcod = ip_int_array("L_IRREPS", moinfo.l_irreps, 3);
  if(errcod == IPE_OK) {
    moinfo.irrep_Rx = moinfo.l_irreps[0];
    moinfo.irrep_Ry = moinfo.l_irreps[1];
    moinfo.irrep_Rz = moinfo.l_irreps[2];
  }
  else {
    fprintf(outfile, "\nYou must supply the irreps of rotation about x, y, and z with the L_IRREPS keyword.\n");
    exit(PSI_RETURN_FAILURE);
  }

  params.maxiter = 50;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-7;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);
  params.diis = 1;
  errcod = ip_boolean("DIIS", &(params.diis), 0);

  if(ip_exist("PROPERTY",0)) {
    errcod = ip_string("PROPERTY", &(params.prop), 0);
    if(strcmp(params.prop,"POLARIZABILITY") && strcmp(params.prop,"ROTATION") && 
       strcmp(params.prop,"ALL")) {
      fprintf(outfile, "Invalid choice of response property: %s\n", params.prop);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.prop = strdup("POLARIZABILITY");

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
    sprintf(local.weakp, "%s", "NONE");
  }

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  if(!strcmp(params.prop,"ALL"))
    fprintf(outfile, "\tProperty        =    POLARIZABILITY + ROTATION\n");
  else
    fprintf(outfile, "\tProperty        =    %s\n", params.prop);
  fprintf(outfile, "\tReference wfn   =    %5s\n",
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tCache Level     =    %1d\n", params.cachelev);
  fprintf(outfile, "\tPrint Level     =    %1d\n",  params.print);
  fprintf(outfile, "\tMaxiter         =    %3d\n",  params.maxiter);
  fprintf(outfile, "\tConvergence     = %3.1e\n", params.convergence);
  fprintf(outfile, "\tDIIS            =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tIrrep X         =    %3s\n", moinfo.labels[moinfo.irrep_x]);
  fprintf(outfile, "\tIrrep Y         =    %3s\n", moinfo.labels[moinfo.irrep_y]);
  fprintf(outfile, "\tIrrep Z         =    %3s\n", moinfo.labels[moinfo.irrep_z]);
  fprintf(outfile, "\tIrrep RX        =    %3s\n", moinfo.labels[moinfo.irrep_Rx]);
  fprintf(outfile, "\tIrrep RY        =    %3s\n", moinfo.labels[moinfo.irrep_Ry]);
  fprintf(outfile, "\tIrrep RZ        =    %3s\n", moinfo.labels[moinfo.irrep_Rz]);
  if(params.omega == 0.0) 
    fprintf(outfile, "\tApplied field   = none\n");
  else 
    fprintf(outfile, "\tApplied field   =    %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	    (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	    _hartree2wavenumbers*params.omega);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
  }
  fprintf(outfile, "\n");
}

