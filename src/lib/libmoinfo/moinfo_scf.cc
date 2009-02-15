#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>

#include "moinfo_scf.h"

extern FILE *outfile;

using namespace std;

namespace psi {

MOInfoSCF::MOInfoSCF() : MOInfoBase()
{
  read_chkpt_data();
  compute_number_of_electrons();
  read_mo_spaces();
  print_mo();
}

MOInfoSCF::~MOInfoSCF() 
{
}

void MOInfoSCF::read_mo_spaces()
{
  /*****************************************************
     See if we're in a subgroup for finite difference
     calculations, by looking to see what OptKing has
     written to the checkpoint file.  Reassign the
     occupation vectors as appropriate.  N.B. the
     SOCC and DOCC are handled by Input (ACS) 
  *****************************************************/

  docc = new int[nirreps];
  actv = new int[nirreps];

  for(int i=0;i<nirreps;i++){
     docc[i] = actv[i] = 0;
  }

  // For single-point geometry optimizations and frequencies
  char *current_displacement_label = chkpt_build_keyword(const_cast<char*>("Current Displacement Irrep"));
  if(chkpt_exist(current_displacement_label)){
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

    int *docc_ref    = new int[nirreps_ref];
    int *actv_ref    = new int[nirreps_ref];

    // build orbital information for current point group
    read_mo_space(nirreps_ref,ndocc,docc_ref,"DOCC");
    read_mo_space(nirreps_ref,nactv,actv_ref,"ACTV ACTIVE SOCC");
    
    for (int h=0; h < nirreps_ref; h++) {
      docc[ correlation[h] ] += docc_ref[h];
      actv[ correlation[h] ] += actv_ref[h];
    }
    
    wfn_sym = correlation[wfn_sym];
    chkpt_set_prefix(save_prefix);
    chkpt_commit_prefix();
    free(save_prefix);
    free(ptgrp_ref);
    delete [] correlation;
    delete [] docc_ref;
    delete [] actv_ref;
  }else{
    // For a single-point only
    read_mo_space(nirreps,ndocc,docc,"DOCC");
    read_mo_space(nirreps,nactv,actv,"ACTV ACTIVE SOCC");
  }

  nactive_ael = nael  - ndocc;
  nactive_bel = nbel  - ndocc;
  
  free(current_displacement_label);
}

void MOInfoSCF::print_mo()
{
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  MOs per irrep:                ");

  for(int i=nirreps;i<8;i++)  
    fprintf(outfile,"     ");
  for(int i=0;i<nirreps;i++)
    fprintf(outfile,"  %s",irr_labs[i]);
  fprintf(outfile," Total");
  fprintf(outfile,"\n  ----------------------------------------------------------------------------");
  print_mo_space(nso,sopi,"Total                         ");
  print_mo_space(ndocc,docc,"Doubly Occupied               ");
  print_mo_space(nactv,actv,"Active/Singly Occupied        ");
  fprintf(outfile,"\n  ----------------------------------------------------------------------------");
  fflush(outfile);
}

}