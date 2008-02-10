/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <iostream>
#include <cmath>
#include "moinfo.h"
#include "utilities.h"
#include "calculation_options.h"

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

MOInfo::MOInfo(int argc, char* argv[],char* id)
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
  
  // The first irrep in the input is 1
  wfn_sym = 1;
  string wavefunction_sym_str = options->get_str_option("WFN_SYM");
  to_lower(wavefunction_sym_str);
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
  root                = options->get_int_option("ROOT") - 1;
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

MOInfo::~MOInfo()
{
  free_memory_info();
  free_memory_mo_spaces();
}

/*!
    \fn MOInfo::read_info()
 */
void MOInfo::read_info()
{
  int value,status;


  /********************
    Get the reference
  ********************/
  char *refstring;
  status = ip_string("REFERENCE",&refstring,0); // Deallocated before next section
  if (status != IPE_OK){
    printf("REFERENCE keyword is missing");
    exit(1);
  }
  else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
    reference = rhf;
  else if (!strcmp(refstring,"UHF"))
    reference = uhf;
  else if (!strcmp(refstring,"ROHF"))
    reference = rohf;
  else if (!strcmp(refstring,"TWOCON"))
    reference = tcscf;
  free(refstring);


  /********************************************************************************
    Get the excitation level, root, print level, and E and T-amp error
  ********************************************************************************/

  // Read the number of electrons or compute it from the CORR_CHARGE parameter
  if(options->get_int_option("NEL")>0){
    nel = options->get_int_option("NEL");
  }else{
    int     natom = chkpt_rd_natom();
    double* zvals = chkpt_rd_zvals();
    nel = 0;
    for(int i=0; i < natom;i++)
      nel += static_cast<int>(zvals[i]);
    nel = nel - options->get_int_option("CORR_CHARGE");
    free(zvals);
  }

  if( ((nel+1-options->get_int_option("MULTP")) % 2) != 0){
    fprintf(outfile,"\n\n\tWrong multiplicity\n\n");
    exit(1);
  }

  nael = (nel + options->get_int_option("MULTP") -1)/2;
  nbel = nel - nael;

  if (false){
    // Don't read these if we're only computing memory requirements
    evals[alpha] = NULL;
    evals[beta]  = NULL;
    clsdpi       = NULL;
    openpi       = NULL;
    orbspi       = chkpt_rd_orbspi();
    sopi         = chkpt_rd_sopi();
  }else{
    // We're doing a real calculation - read the checkpoint file
    /***********************************
      Read the Fock matrix eigenvalues
    ***********************************/
    switch(reference){
    case rhf:
      evals[alpha]  = chkpt_rd_evals();
      evals[beta]   = chkpt_rd_evals();
      break; //end of rhf case
    case uhf:
      evals[alpha]  = chkpt_rd_alpha_evals();
      evals[beta]   = chkpt_rd_beta_evals();
      break; //end of uhf case
    case rohf:
      evals[alpha]  = chkpt_rd_evals();
      evals[beta]   = chkpt_rd_evals();
      break; //end of rohf case
    case tcscf:
      evals[alpha]  = chkpt_rd_evals();
      evals[beta]   = chkpt_rd_evals();
      break; //end of rohf case
    default:
      printf("REFERENCE %s not implemented yet in the MOInfo class",refstring);
      exit(1);
    }
    clsdpi         = chkpt_rd_clsdpi();
    openpi         = chkpt_rd_openpi();
    scf_energy     = chkpt_rd_escf();
    orbspi         = chkpt_rd_orbspi();
    sopi           = chkpt_rd_sopi();
  }
  /***********************************
    Read Nuclear,SCF and other stuff
  ***********************************/
  nuclear_energy = chkpt_rd_enuc();
  nso            = chkpt_rd_nso();
  norbs          = chkpt_rd_nmo();
  nirreps        = chkpt_rd_nirreps();
  irr_labs       = chkpt_rd_irr_labs();
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

  for(int i=0;i<nirreps;i++){
     focc[i]=docc[i]=actv[i]=fvir[i]=0;
  }

  // For single-point geometry optimizations and frequencies
  if(chkpt_exist(chkpt_build_keyword("Current Displacement Irrep"))){
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


void MOInfo::read_mo_space(int nirreps_ref,int& n, int* mo, char* label)
{
  int size;
  int stat=ip_count(label,&size,0);
  n=0;
  if(stat == IPE_OK){
    if(size==nirreps_ref){
      for(int i=0;i<size;i++){
        stat=ip_data(label,"%d",(&(mo)[i]),1,i);
        n+=mo[i];
      }
    }else{
      fprintf(outfile,"\n\nThe size of the %s array (%d) does not match the number of irreps (%d), please fix the input file",label,size,nirreps_ref);
      fflush(outfile);
      exit(1);
    }
  }else{ // The keyword label was not found
    for(int i=0;i<nirreps_ref;i++)
      mo[i]=0;
    n = 0;
  }
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
  print_mo_space(navir,avir,"Active Virtual                  ");
  print_mo_space(nfvir,fvir,"Frozen Virtual                  ");
  fflush(outfile);
}


void MOInfo::print_mo_space(int& n, int* mo, char* label)
{
  fprintf(outfile,"\n  %s",label);

  for(int i=nirreps;i<8;i++)
    fprintf(outfile,"     ");
  for(int i=0;i<nirreps;i++)
    fprintf(outfile," %3d ",mo[i]);
  fprintf(outfile,"  %3d",n);
}


void MOInfo::write_scf_mos()
{
  chkpt_wt_scf(scf);
}

void MOInfo::correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& arr)
{ /* This is a hack from input!!! (ACS) */
  int  i;

  if (strcmp(ptgrp,"C1 ") == 0)
    nirreps_old = 1;
  else if (strcmp(ptgrp,"Cs ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"Ci ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2 ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2v") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2 ") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"C2h") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2h") == 0)
    nirreps_old = 8;
  else {
    fprintf(outfile,"point group %s unknown.\n",ptgrp);
    exit(1);
  }

  arr = new int[nirreps_old];

  if (irrep == 0) { /* return identity */
    nirreps_new = nirreps_old;
    for (i=0; i<nirreps_old; ++i)
      arr[i] = i;
    return;
  }

  nirreps_new = nirreps_old / 2;
  if ((strcmp(ptgrp,"C1 ") == 0) || (strcmp(ptgrp,"Cs ") == 0) ||
      (strcmp(ptgrp,"Ci ") == 0) || (strcmp(ptgrp,"C2 ") == 0) ) {
        arr[0] = 0; arr[1] = 0;
  }
  else if ( (strcmp(ptgrp,"C2v") == 0) || (strcmp(ptgrp,"D2 ") == 0) ||
            (strcmp(ptgrp,"C2h") == 0) ) {
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
    }
  }
  else if (strcmp(ptgrp,"D2h") == 0) {
    /* 1,2,3 give C2h displaced geometries */
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 2; arr[6] = 3;  arr[7] = 3;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 3; arr[6] = 2;  arr[7] = 3;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
      arr[4] = 2;  arr[5] = 3; arr[6] = 3;  arr[7] = 2;
    }
    /* 4 gives D2 displaced geometries */
    else if (irrep == 4) { /* D2 */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 0;  arr[5] = 1; arr[6] = 2;  arr[7] = 3;
    }
    /* displacements along irreps 5,6,7 make C2v structures */
    /* care is taken to make sure definition of b1 and b2 will
       match those that input will generate - the following seems to work:
       b1u disp: has C2(z), b2 irrep symmetric wrt sigma(yz)
       b2u disp: has C2(y), b2 irrep symmetric wrt sigma(xy)
       b3u disp: has C2(x), b2 irrep symmetric wrt sigma(xz) */
    else if (irrep == 5) { /* b1u */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 1;  arr[5] = 0; arr[6] = 3;  arr[7] = 2;
    }
    else if (irrep == 6) { /* b2u */
      arr[0] = 0;  arr[1] = 3; arr[2] = 1;  arr[3] = 2;
      arr[4] = 1;  arr[5] = 2; arr[6] = 0;  arr[7] = 3;
    }
    else if (irrep == 7) { /* b3u */
      arr[0] = 0;  arr[1] = 2; arr[2] = 3;  arr[3] = 1;
      arr[4] = 1;  arr[5] = 3; arr[6] = 2;  arr[7] = 0;
    }
  }
  else {
    fprintf(outfile,"Point group unknown for correlation table.\n");
  }

  return;
}

/*!
    \fn MOInfo::free_memory_info()
 */
void MOInfo::free_memory_info()
{
  /// @todo implement me
  for(int i=0;i<nirreps;i++)
     free(irr_labs[i]);
  free(irr_labs);
  if(evals[alpha]!=NULL)
    free(evals[alpha]);
  if(evals[beta]!=NULL)
    free(evals[beta]);
}

/*!
    \fn MOInfo::free_memory_mo_spaces()
 */
void MOInfo::free_memory_mo_spaces()
{
  /// @todo implement me
  free(orbspi);
  free(sopi);
  if(clsdpi != NULL);
    free(clsdpi);
  if(openpi != NULL);
    free(openpi);
  if(scf != NULL);
    free_block(scf);
  for(int i=0;i<nirreps;i++)
    free_block(scf_irrep[i]);
  delete[] scf_irrep;
  delete[] ioff;
  delete[] focc;
  delete[] docc;
  delete[] actv;
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

}} /* End Namespaces */