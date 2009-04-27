// Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

// STL
#include <numeric>

// PSI Libraries
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libutil/libutil.h>
#include <libqt/qt.h>

#include "moinfo.h"

extern FILE *outfile;

using namespace std;

namespace psi {

MOInfo::MOInfo() : MOInfoBase()
{
  /***************
    Set defaults
  ***************/
  no_damp_convergence = 1.0e-9;
  dgemm_timing        = 0.0;
  scf                 = NULL;

  nfocc = 0;
  nfvir = 0;
  nactv_docc = 0;
  nocc  = 0;
  nvir  = 0;
  nall  = 0;
  nextr = 0;

  read_info();
  read_mo_spaces();
  compute_mo_mappings();
  print_info();
  print_mo();

  // Determine the wave function irrep

  // The first irrep in the input is 1
  wfn_sym = 1;
  string wavefunction_sym_str = options_get_str("WFN_SYM");
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
  root = options_get_int("ROOT") - 1;
}

MOInfo::~MOInfo()
{
  free_memory();
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


void MOInfo::read_info()
{
  /***********************************
    Read Nuclear,SCF and other stuff
  ***********************************/
  read_chkpt_data();
  nmo            = _default_chkpt_lib_->rd_nmo();
  compute_number_of_electrons();
  scf_energy     = _default_chkpt_lib_->rd_escf();
  mopi           = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_orbspi());
  scf            = _default_chkpt_lib_->rd_scf();
  scf_irrep      = new double**[nirreps];
  for(int i=0;i<nirreps;i++)
    scf_irrep[i] = _default_chkpt_lib_->rd_scf_irrep(i);
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

  focc.assign(nirreps,0);
  docc.assign(nirreps,0);
  actv.assign(nirreps,0);
  fvir.assign(nirreps,0);
  extr.assign(nirreps,0);
  occ.assign(nirreps,0);
  vir.assign(nirreps,0);
  all.assign(nirreps,0);
  actv_docc.assign(nirreps,0);

  // For single-point geometry optimizations and frequencies
  if(_default_chkpt_lib_->exist(_default_chkpt_lib_->build_keyword(const_cast<char *>("Current Displacement Irrep")))){
    int   disp_irrep  = _default_chkpt_lib_->rd_disp_irrep();
    char *save_prefix = _default_chkpt_lib_->rd_prefix();
    int nirreps_ref;

    // read symmetry info and MOs for undisplaced geometry from
    // root section of checkpoint file
    _default_chkpt_lib_->reset_prefix();
    _default_chkpt_lib_->commit_prefix();

    char *ptgrp_ref = _default_chkpt_lib_->rd_sym_label();

    // Lookup irrep correlation table
    int* correlation;
    correlate(ptgrp_ref, disp_irrep, nirreps_ref, nirreps,correlation);

    intvec focc_ref;
    intvec docc_ref;
    intvec actv_ref;
    intvec fvir_ref;
    intvec actv_docc_ref;

    // build orbital information for current point group
    // Read the values stored in the chkpt
    // override if the user defines values
    focc_ref = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_frzcpi());
    docc_ref = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_clsdpi());
    actv_ref = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_openpi());

    for (int h = 0; h < nirreps_ref; h++)
      docc_ref[h] -= focc_ref[h];

    nfocc = std::accumulate( focc_ref.begin(), focc_ref.end(), 0 );
    ndocc = std::accumulate( docc_ref.begin(), docc_ref.end(), 0 );
    nactv = std::accumulate( actv_ref.begin(), actv_ref.end(), 0 );

    read_mo_space(nirreps_ref,nfocc,focc_ref,"CORR_FOCC FROZEN_DOCC",false);
    read_mo_space(nirreps_ref,ndocc,docc_ref,"CORR_DOCC RESTRICTED_DOCC",false);
    read_mo_space(nirreps_ref,nactv,actv_ref,"CORR_ACTV ACTIVE ACTV",false);

    read_mo_space(nirreps_ref,nfvir,fvir_ref,"CORR_FVIR FROZEN_UOCC");
    read_mo_space(nirreps_ref,nactv_docc,actv_docc_ref,"ACTIVE_DOCC");

    for (int h = 0; h < nirreps_ref; h++) {
      focc[ correlation[h] ]      += focc_ref[h];
      docc[ correlation[h] ]      += docc_ref[h];
      actv[ correlation[h] ]      += actv_ref[h];
      fvir[ correlation[h] ]      += fvir_ref[h];
      actv_docc[ correlation[h] ] += actv_docc_ref[h];
    }
    wfn_sym = correlation[wfn_sym];
    _default_chkpt_lib_->set_prefix(save_prefix);
    _default_chkpt_lib_->commit_prefix();
    free(save_prefix);
    free(ptgrp_ref);
    delete [] correlation;
  }else{
    // For a single-point only
    focc = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_frzcpi());
    docc = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_clsdpi());
    actv = read_chkpt_intvec(nirreps,_default_chkpt_lib_->rd_openpi());

    for (int h = 0; h < nirreps; h++)
      docc[h] -= focc[h];

    nfocc = std::accumulate( focc.begin(), focc.end(), 0 );
    ndocc = std::accumulate( docc.begin(), docc.end(), 0 );
    nactv = std::accumulate( actv.begin(), actv.end(), 0 );

    read_mo_space(nirreps,nfocc,focc,"CORR_FOCC FROZEN_DOCC",false);
    read_mo_space(nirreps,ndocc,docc,"CORR_DOCC RESTRICTED_DOCC",false);
    read_mo_space(nirreps,nactv,actv,"CORR_ACTV ACTIVE ACTV",false);
    read_mo_space(nirreps,nfvir,fvir,"CORR_FVIR FROZEN_UOCC");
    read_mo_space(nirreps,nactv_docc,actv_docc,"ACTIVE_DOCC");
  }

  // Compute the number of external orbitals
  nextr = 0;
  for(int h = 0; h < nirreps; ++h){
     extr[h]= mopi[h] - focc[h] - docc[h] - actv[h] - fvir[h];
     occ[h] = docc[h] + actv[h];
     vir[h] = actv[h] + extr[h];
     all[h] = mopi[h] - focc[h] - fvir[h];
     nextr += extr[h];
  }
  nall        = nmo  - nfocc - nfvir;
  nactive_ael = nael - ndocc - nfocc;
  nactive_bel = nbel - ndocc - nfocc;
  nocc        = ndocc + nactv;
  nvir        = nactv + nextr;

  /*********************************************
    Define the symmetry of each  non-frozen MO
  **********************************************/
  all_sym.resize(nall);
  int index_mo  = 0;
  for(int h = 0; h < nirreps; ++h){
     for(int i  =0; i < all[i]; ++i){
      all_sym[index_mo] = h;
      index_mo++;
    }
  }

  /***************************************************************
    Build the array that connects the non-frozen MOs (all) to the
    the complete list of MOs (mo). Used when frozen MOs are used.
  ****************************************************************/
  all_to_mo.resize(nall);
  int index_all = 0;
  index_mo  = 0;
  for(int h = 0; h < nirreps; ++h){
    index_mo += focc[h];
    for(int i = 0; i < all[h]; ++i){
      all_to_mo[index_all]=index_mo;
      index_all++;
      index_mo++;
    }
    index_mo += fvir[h];
  }

  // The mapping of the MOs (mo) to the non-frozen MOs (all)
  // Set size to nmo and all elements to -1
  mo_to_all.assign(nmo,-1);
  // Set the mappings
  for(int i = 0; i < nall; ++i)
    mo_to_all[all_to_mo[i]]=i;
}

/**
    MOInfo::print_mo_spaces()
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
  print_mo_space(nmo,mopi,"Total                           ");
  print_mo_space(nfocc,focc,"Frozen Occupied                 ");
  print_mo_space(ndocc,docc,"Doubly Occupied                 ");
  print_mo_space(nactv,actv,"Active                          ");
  if(nactv_docc > 0){
    print_mo_space(nactv_docc,actv_docc,"Active Doubly Occupied          ");
  }
  print_mo_space(nextr,extr,"External                        ");
  print_mo_space(nfvir,fvir,"Frozen Virtual                  ");
  fflush(outfile);
}

/**
 *   MOInfo::free_memory()
 */
void MOInfo::free_memory()
{
  if(scf != NULL);
    free_block(scf);
  for(int i=0;i<nirreps;i++)
    free_block(scf_irrep[i]);
  delete[] scf_irrep;
}

}
