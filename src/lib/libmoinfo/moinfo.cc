// Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

// PSI Libraries
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libutil/libutil.h>
#include <libqt/qt.h>

#include "moinfo.h"

extern FILE *outfile;

using namespace std;

MOInfo::MOInfo() : MOInfoBase()
{
  /***************
    Set defaults
  ***************/
  no_damp_convergence = 1.0e-9;
  dgemm_timing        = 0.0;
  scf                 = NULL;
//   tuning();
  read_info();
  read_mo_spaces();
  compute_mo_mappings();
  print_info();
  print_mo();
  
  // The first irrep in the input is 1 // TODO: Add control on the wfn_sym keyword
  wfn_sym = 1;
  string wavefunction_sym_str = options_get_str("WFN_SYM");
  to_lower(wavefunction_sym_str);
  
//  vector<string> irr_labels_vector;
//  for(int h = 0; h < nirreps; ++h){
//    string irr_label_str = irr_labs[h];
//    irr_labels_vector.push_back(to_lower(irr_label_str));
//  }
  
  for(int h=0; h < nirreps; ++h){
    string irr_label_str = irr_labs[h];
    trim_spaces(irr_label_str);
    to_lower(irr_label_str);
    if(wavefunction_sym_str==irr_label_str){
      wfn_sym = h;
    }
    if(to_string(h+1) == wavefunction_sym_str){
      wfn_sym = h;
    }
  }

  // The lowest root in the input is 1, here we subtract one
  root                = options_get_int("ROOT") - 1;
}

MOInfo::~MOInfo()
{
  free_memory_info();
  free_memory_mo_spaces();
}

void MOInfo::setup_model_space()
{
  // NB This code could be places elsewhere
  references.clear();
  alpha_internal_excitations.clear();
  beta_internal_excitations.clear();
  sign_internal_excitations.clear();
  all_refs.clear();
  unique_refs.clear();
  closed_shell_refs.clear();
  unique_open_shell_refs.clear();

  build_model_space();
  print_model_space();
  make_internal_excitations();
}


/*!
    \fn MOInfo::read_info()
 */
void MOInfo::read_info()
{
  int value,status;


//  /********************
//    Get the reference
//  ********************/
//  char *refstring;
//  status = ip_string("REFERENCE",&refstring,0); // Deallocated before next section
//  if (status != IPE_OK){
//    printf("REFERENCE keyword is missing");
//    exit(1);
//  }
//  else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
//    reference = rhf;
//  else if (!strcmp(refstring,"UHF"))
//    reference = uhf;
//  else if (!strcmp(refstring,"ROHF"))
//    reference = rohf;
//  else if (!strcmp(refstring,"TWOCON"))
//    reference = tcscf;
//  free(refstring);


  /********************************************************************************
    Get the excitation level, root, print level, and E and T-amp error
  ********************************************************************************/

  // Read the number of electrons or compute it from the CORR_CHARGE parameter

  if (false){
    // Don't read these if we're only computing memory requirements
    evals[alpha] = NULL;
    evals[beta]  = NULL;
    clsdpi       = NULL;
    openpi       = NULL;
    orbspi       = chkpt_rd_orbspi();
//    sopi         = chkpt_rd_sopi();
  }else{
    // We're doing a real calculation - read the checkpoint file
    /***********************************
      Read the Fock matrix eigenvalues
    ***********************************/
//    switch(reference){
//    case rhf:
//      evals[alpha]  = chkpt_rd_evals();
//      evals[beta]   = chkpt_rd_evals();
//      break; //end of rhf case
//    case uhf:
//      evals[alpha]  = chkpt_rd_alpha_evals();
//      evals[beta]   = chkpt_rd_beta_evals();
//      break; //end of uhf case
//    case rohf:
//      evals[alpha]  = chkpt_rd_evals();
//      evals[beta]   = chkpt_rd_evals();
//      break; //end of rohf case
//    case tcscf:
//      evals[alpha]  = chkpt_rd_evals();
//      evals[beta]   = chkpt_rd_evals();
//      break; //end of rohf case
//    default:
//      printf("REFERENCE %s not implemented yet in the MOInfo class",refstring);
//      exit(1);
//    }
    clsdpi         = chkpt_rd_clsdpi();
    openpi         = chkpt_rd_openpi();
    scf_energy     = chkpt_rd_escf();
    orbspi         = chkpt_rd_orbspi();
  }
  /***********************************
    Read Nuclear,SCF and other stuff
  ***********************************/
  
  read_chkpt_data();
  compute_number_of_electrons();

  norbs          = chkpt_rd_nmo();
  scf            = chkpt_rd_scf();
  scf_irrep      = new double**[nirreps];
  for(int i=0;i<nirreps;i++)
    scf_irrep[i]  = chkpt_rd_scf_irrep(i);
}


/*!
    \fn MOInfo::print_info()
 */
void MOInfo::print_info()
{
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  ==============================================================================");
  fprintf(outfile,"\n  System Info:");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  fprintf(outfile,"\n  Nuclear Energy   = %-15.9f  SCF Energy       = %-15.9f",nuclear_energy,scf_energy);
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  MOs and Symmetry:");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  fprintf(outfile,"\n  nirreps          = %-10d",nirreps);
  fprintf(outfile,"\n  nso              = %-10d       nmo              = %-10d",nso,nmo);
  fprintf(outfile,"\n  nael             = %-10d       nbel             = %-10d",nael,nbel);
  fprintf(outfile,"\n  nactive_ael      = %-10d       nactive_bel      = %-10d",nactive_ael,nactive_bel);
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  Details of the Computation:");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
}

/*!
    \fn MOInfo::read_mo_spaces()
 */
void MOInfo::read_mo_spaces()
{
  /*****************************************************
     See if we're in a subgroup for finite difference
     calculations, by looking to see what OptKing has
     written to the checkpoint file.  Reassign the
     occupation vectors as appropriate.  N.B. the
     SOCC and DOCC are handled by Input (ACS) 
  *****************************************************/

  focc = new int[nirreps];
  docc = new int[nirreps];
  actv = new int[nirreps];
  fvir = new int[nirreps];
  avir = new int[nirreps];
  occ  = new int[nirreps];
  vir  = new int[nirreps];
  actv_docc = NULL;

  for(int i=0;i<nirreps;i++){
     focc[i]=docc[i]=actv[i]=fvir[i]=0;
  }

  // For single-point geometry optimizations and frequencies
  if(chkpt_exist(chkpt_build_keyword(const_cast<char *>("Current Displacement Irrep")))){
    int   disp_irrep  = chkpt_rd_disp_irrep();
    char *save_prefix = chkpt_rd_prefix();
    int nirreps_ref;

    // read symmetry info and MOs for undisplaced geometry from
    // root section of checkpoint file
    chkpt_reset_prefix();
    chkpt_commit_prefix();

    char *ptgrp_ref = chkpt_rd_sym_label();

    // Lookup irrep correlation table
    int* correlation;
    correlate(ptgrp_ref, disp_irrep, nirreps_ref, nirreps,correlation);

    int *focc_ref    = new int[nirreps_ref];
    int *docc_ref    = new int[nirreps_ref];
    int *actv_ref    = new int[nirreps_ref];
    int *fvir_ref    = new int[nirreps_ref];

    // build orbital information for current point group
    read_mo_space(nirreps_ref,nfocc,focc_ref,"CORR_FOCC");
    read_mo_space(nirreps_ref,ndocc,docc_ref,"CORR_DOCC");
    read_mo_space(nirreps_ref,nactv,actv_ref,"CORR_ACTV");
    read_mo_space(nirreps_ref,nfvir,fvir_ref,"CORR_FVIR");
     
    
    for (int h=0; h < nirreps_ref; h++) {
      focc[ correlation[h] ] += focc_ref[h];
      docc[ correlation[h] ] += docc_ref[h];
      actv[ correlation[h] ] += actv_ref[h];
      fvir[ correlation[h] ] += fvir_ref[h];
    }
    wfn_sym = correlation[wfn_sym];
    chkpt_set_prefix(save_prefix);
    chkpt_commit_prefix();
    free(save_prefix);
    free(ptgrp_ref);
    delete [] correlation;
    delete [] focc_ref;
    delete [] docc_ref;
    delete [] actv_ref;
    delete [] fvir_ref;
  }else{
    // For a single-point only
    read_mo_space(nirreps,nfocc,focc,"CORR_FOCC");
    read_mo_space(nirreps,ndocc,docc,"CORR_DOCC");
    read_mo_space(nirreps,nactv,actv,"CORR_ACTV");
    read_mo_space(nirreps,nfvir,fvir,"CORR_FVIR");
    if(options_get_str("CORR_WFN") == "CCSD_MP2"){
      actv_docc = new int[nirreps];
      read_mo_space(nirreps,nactv_docc,actv_docc,"CORR_ACTV_DOCC");
    }
  }

  // Compute the number of active virtuals
  navir=0;
  for(int i=0;i<nirreps;i++){
     avir[i]=orbspi[i]-focc[i]-docc[i]-actv[i]-fvir[i];
     occ[i]=docc[i]+actv[i];
     vir[i]=actv[i]+avir[i];
     navir+=avir[i];
  }
  nmo         = norbs - nfocc - nfvir;
  nactive_ael = nael  - ndocc - nfocc;
  nactive_bel = nbel  - ndocc - nfocc;
  nocc        = ndocc + nactv;
  nvir        = nactv + navir;

  /*********************************
    Define the symmetry of each MO
  *********************************/
  mo_irr = new int[nmo];
  int mo_index=0;
  for(int i=0;i<nirreps;i++){
     for(int j=0;j<orbspi[i]-focc[i]-fvir[i];j++)
      mo_irr[mo_index++]=i;
  }

  /********************************************************
    Build the array that connects the non-frozen MOs to the
    the complete list of MOs. Used when frozen MOs are used.
  ********************************************************/
  nonfrozen_to_all = new int[nmo];
  int index_all =0;
  int index_nonfrozen =0;
  for(int irrep=0;irrep<nirreps;irrep++){
    int nonfrozen = orbspi[irrep]-focc[irrep]-fvir[irrep];
    index_all+=focc[irrep];
    for(int i=0;i<nonfrozen;i++)
      nonfrozen_to_all[index_nonfrozen++]=index_all++;
    index_all+=fvir[irrep];
  }

  all_to_nonfrozen = new int[norbs];
  for(int i=0;i<norbs;i++)
    all_to_nonfrozen[i]=-1;
  for(int i=0;i<nmo;i++)
    all_to_nonfrozen[nonfrozen_to_all[i]]=i;
}

/*!
    \fn MOInfo::print_mo_spaces()
 */
void MOInfo::print_mo()
{
  /// @todo implement me
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  MOs per irrep:                  ");

  for(int i=nirreps;i<8;i++)  
    fprintf(outfile,"     ");
  for(int i=0;i<nirreps;i++)
    fprintf(outfile,"  %s",irr_labs[i]);
  fprintf(outfile," Total");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  print_mo_space(nso,orbspi,"Total                           ");
  print_mo_space(nfocc,focc,"Frozen Occupied                 ");
  print_mo_space(ndocc,docc,"Doubly Occupied                 ");
  print_mo_space(nactv,actv,"Active                          ");

  if(options_get_str("CORR_WFN") == "CCSD_MP2"){  
    print_mo_space(nactv_docc,actv_docc,"Active Doubly Occupied          ");
  }

  print_mo_space(navir,avir,"Active Virtual                  ");
  print_mo_space(nfvir,fvir,"Frozen Virtual                  ");
  fflush(outfile);
}



/*!
    \fn MOInfo::free_memory_info()
 */
void MOInfo::free_memory_info()
{
  /// @todo implement me
//  if(evals[alpha]!=NULL)
//    free(evals[alpha]);
//  if(evals[beta]!=NULL)
//    free(evals[beta]);
}

/*!
    \fn MOInfo::free_memory_mo_spaces()
 */
void MOInfo::free_memory_mo_spaces()
{
  /// @todo implement me
  free(orbspi);
  if(clsdpi != NULL);
    free(clsdpi);
  if(openpi != NULL);
    free(openpi);
  if(scf != NULL);
    free_block(scf);
  for(int i=0;i<nirreps;i++)
    free_block(scf_irrep[i]);
  if(actv_docc != NULL)
    delete[] actv_docc;
  delete[] scf_irrep;
  delete[] focc;
  delete[] avir;
  delete[] fvir;
  delete[] occ;
  delete[] vir;
  delete[] mo_irr;
  delete[] nonfrozen_to_all;
  delete[] all_to_nonfrozen;
  delete[] so_to_pitzer;
  delete[] orbs_to_pitzer;
  delete[] docc_to_pitzer;
  delete[] act_to_pitzer;
  delete[] ext_to_pitzer;
  delete[] occ_to_pitzer;
  delete[] vir_to_pitzer;
  delete[] all_to_pitzer;
  delete[] occ_to_vir;
  delete[] all_to_occ;
  delete[] all_to_vir;
  delete[] act_to_occ;
  delete[] act_to_vir;
  delete[] occ_to_act;
  delete[] vir_to_act;
  delete[] first_orbs_pitzer;
  delete[] last_orbs_pitzer;
  delete[] first_so_pitzer;
  delete[] last_so_pitzer;
  delete[] first_occupied_pitzer[alpha];
  delete[] last_occupied_pitzer[alpha];
  delete[] first_virtual_pitzer[alpha];
  delete[] last_virtual_pitzer[alpha];
  delete[] first_active_pitzer[alpha];
  delete[] last_active_pitzer[alpha];

  delete[] is_act_in_occ;
  delete[] is_act_in_vir;
}
