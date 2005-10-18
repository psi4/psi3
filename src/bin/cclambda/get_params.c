#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void get_params(void)
{
  int errcod, iconv,i,j,k,l,prop_sym,prop_root;
  char lbl[32];
  char *junk;

  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "CCSD") && strcmp(params.wfn, "CCSD_T") &&
     strcmp(params.wfn, "EOM_CCSD") && strcmp(params.wfn, "LEOM_CCSD") &&
     strcmp(params.wfn, "BCCD") && strcmp(params.wfn,"BCCD_T") &&
     strcmp(params.wfn, "CC2") && strcmp(params.wfn,"CC3") &&
     strcmp(params.wfn, "EOM_CC2") && strcmp(params.wfn, "EOM_CC3")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }

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

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

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

  params.diis = 1;
  errcod = ip_boolean("DIIS", &params.diis, 0);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);
  params.aobasis = 0;  /* AO basis code not yet working for lambda */

  /* Here is the current overall logic:
     if --all, do all states included ground state
     else if (params.ground), do just ground state
     else if --zeta {
     L_irr from CC_INFO (A1 for gradients)
     labels for zeta
     }
     else {
     if prop_sym is in input use prop_sym and prop_root
     else compute L for last state requested in input
     }
  */

  /* count number of states to converge */
  if (params.all) { /* compute LAMPS for all Rs plus ground state */
    if (ip_exist("STATES_PER_IRREP",0)) {
      ip_count("STATES_PER_IRREP", &i, 0);
      if (i != moinfo.nirreps) {
        fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps) ;
        exit(0);
      }
      for (i=0;i<moinfo.nirreps;++i) {
        ip_data("STATES_PER_IRREP","%d",&j,1,i);
        params.nstates += j;
      }
    }
    params.nstates += 1; /* for ground state */
  }
  else {
    params.nstates = 1;
  }

  pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));

  if (params.all) {
    /* for ground state */
    pL_params[0].irrep = 0;
    pL_params[0].R0 = 1.0;
    pL_params[0].cceom_energy = 0.0;
    pL_params[0].ground = 1;
    pL_params[0].root = -1;
    sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
    sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
    sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
    sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
    sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
    sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);

    l=1; /* index of root */
    if (ip_exist("STATES_PER_IRREP",0)) {
      for (i=0;i<moinfo.nirreps;++i) {
        ip_data("STATES_PER_IRREP","%d",&j,1,i);
        for (k=0;k<j;++k) {
          pL_params[l].irrep = i^moinfo.sym;

          if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
            sprintf(lbl,"EOM CC2 Energy for root %d %d", pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].cceom_energy),sizeof(double));
            sprintf(lbl,"EOM CC2 R0 for root %d %d",pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].R0),sizeof(double));
          }
          else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
            sprintf(lbl,"EOM CCSD Energy for root %d %d", pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].cceom_energy),sizeof(double));
            sprintf(lbl,"EOM CCSD R0 for root %d %d",pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].R0),sizeof(double));
          }
          else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
            sprintf(lbl,"EOM CC3 Energy for root %d %d", pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].cceom_energy),sizeof(double));
            sprintf(lbl,"EOM CC3 R0 for root %d %d",pL_params[l].irrep, k);
            psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[l].R0),sizeof(double));
          }

          sprintf(pL_params[l].L1A_lbl,"LIA %d %d",pL_params[l].irrep, k);
          sprintf(pL_params[l].L1B_lbl,"Lia %d %d",pL_params[l].irrep, k);
          sprintf(pL_params[l].L2AA_lbl,"LIJAB %d %d",pL_params[l].irrep, k);
          sprintf(pL_params[l].L2BB_lbl,"Lijab %d %d",pL_params[l].irrep, k);
          sprintf(pL_params[l].L2AB_lbl,"LIjAb %d %d",pL_params[l].irrep, k);
          sprintf(pL_params[l].L2RHF_lbl,"2LIjAb - LIjbA %d %d",pL_params[l].irrep, k);
          pL_params[l].root = k;
          pL_params[l++].ground = 0;
        }
      }
    }
  }
  else if (params.zeta) {
    psio_read_entry(CC_INFO, "XI Irrep", (char *) &i,sizeof(int));
    fprintf(outfile,"\tIrrep of Zeta       (CC_INFO) = %d\n", i);
    prop_root = 0;
    prop_sym = i;
    pL_params[0].irrep = prop_sym;
    pL_params[0].root = prop_root;
    pL_params[0].ground = 0;
    pL_params[0].cceom_energy = 0.0; /* don't want energy in denominator */
    pL_params[0].R0 = 0.0; /* <Zeta0|R0> = 0, since zeta_0 = 0 */
    sprintf(pL_params[0].L1A_lbl,"ZIA");
    sprintf(pL_params[0].L1B_lbl,"Zia");
    sprintf(pL_params[0].L2AA_lbl,"ZIJAB");
    sprintf(pL_params[0].L2BB_lbl,"Zijab");
    sprintf(pL_params[0].L2AB_lbl,"ZIjAb");
    sprintf(pL_params[0].L2RHF_lbl,"2ZIjAb - ZIjbA");
  }
  else if (params.ground) {
    pL_params[0].irrep = 0;
    pL_params[0].R0 = 1.0;
    pL_params[0].cceom_energy = 0.0;
    pL_params[0].ground = 1;
    sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
    sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
    sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
    sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
    sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
    sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);
    pL_params[0].root = -1;
  }
  else { /* use input to determine which L to get */
    if (ip_exist("PROP_SYM",0)) {
      ip_data("PROP_SYM","%d",&(prop_sym),0);
      prop_sym -= 1;
      prop_sym = moinfo.sym^prop_sym;
      ip_data("PROP_ROOT","%d",&(prop_root),0);
      prop_root -= 1;
      pL_params[0].irrep = prop_sym;
      pL_params[0].ground = 0;

      if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
        sprintf(lbl,"EOM CC3 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC3 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }

      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",prop_sym, prop_root);
      pL_params[0].root = prop_root;
    }
    else { /* use last root requested */
      for (i=0;i<moinfo.nirreps;++i) {
        j=0;
        ip_data("STATES_PER_IRREP","%d",&j,1,i);
        if (j>0) {
          prop_root = j-1;
          prop_sym = i^moinfo.sym;
        }
      }
      pL_params[0].irrep = prop_sym;
      pL_params[0].ground = 0;
      sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
      psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
      sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
      psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",prop_sym, prop_root);
      pL_params[0].root = prop_root;
    }
  }

  if(ip_exist("ABCD",0)) {
    errcod = ip_string("ABCD", &(params.abcd), 0);
    if(strcmp(params.abcd,"NEW") && strcmp(params.abcd,"OLD")) {
      fprintf(outfile, "Invalid ABCD algorithm: %s\n", params.abcd);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.abcd = strdup("NEW");

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
  
  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;
  ip_boolean("LOCAL_FILTER_SINGLES", &(local.filter_singles), 0);

  local.cphf_cutoff = 0.10;
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
  else if(params.local && params.dertype == 3)
    local.pairdef = strdup("RESPONSE");
  else if(params.local)
    local.pairdef = strdup("BP");

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter       =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence   = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart       =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\tCache Level   =     %1d\n", params.cachelev);
  fprintf(outfile, "\tDIIS          =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tAO Basis      =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tABCD            =     %s\n", params.abcd);
  fprintf(outfile, "\tExcited State  =     %s\n", 
          params.ground ? "No" : "Yes");
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef);
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }

  fprintf(outfile,"\tParamaters for left-handed eigenvectors:\n");
  fprintf(outfile,"\tIrr   Root  Ground-State?    EOM energy        R0\n");
  for (i=0; i<params.nstates; ++i) {
    fprintf(outfile,"\t%3d %5d %13s %14.10lf %14.10lf\n", pL_params[i].irrep, pL_params[i].root,
	    (pL_params[i].ground ? "Yes":"No"), pL_params[i].cceom_energy, pL_params[i].R0);
  }

  /*
  for (i=0; i<params.nstates; ++i) {
    fprintf(outfile,"\tLabels for state %d:\n\t%s, %s, %s, %s, %s, %s\n",
	    i,pL_params[i].L1A_lbl,pL_params[i].L1B_lbl,
            pL_params[i].L2AA_lbl,pL_params[i].L2BB_lbl,
	    pL_params[i].L2AB_lbl, pL_params[i].L2RHF_lbl);
  }
  */

  return;
}
