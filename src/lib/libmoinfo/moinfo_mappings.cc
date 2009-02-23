#include "moinfo.h"

extern FILE *outfile;

using namespace std;

namespace psi {

/*!
 * \fn MOInfo::compute_mo_mappings()
 */
void MOInfo::compute_mo_mappings()
{
//
//  ------------------------------------------------------------------------------
//  |0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15|16 17 18 19 20 21 22 23 24 25 26 27 28|
//  ------------------------------------------------------------------------------
//  |               Irrep 0               |               Irrep 1                |
//  | focc  | docc  |act|    extr   |fvir |focc | docc   | act |    extr   | fvir|
//  |-------|     occ   |-----------------------|       occ    |-----------------|
//  |---------------|     vir       |--------------------|     vir         |-----|
//  |-------|          all          |-----------|            all           |-----|
//  |                                     mo                                     |
//  |----------------------------------------------------------------------------|
//

  /********************************************************
    Build the arrays that connect the subspaces of the
    non-frozen MOs to the the non-frozen list of MOs
    in Pitzer. This is what you need in all computations.

    For the example above:

    docc_to_mo = (4 5 6 7 18 19 20)
    actv_to_mo = (8 9 21 22)
    extr_to_mo = (10 11 12 13 23 24 25 26)
    occ_to_mo  = (4 5 6 7 8 9 18 19 20 21 22)
    vir_to_mo  = (8 9 10 11 12 13 21 22 23 24 25 26)
    all_to_mo  = (4 5 6 7 8 9 10 11 12 13 18 19 20
                      21 22 23 24 25 26)
    all_to_occ     = (0 1 2 3 4 5 -1 -1 -1 -1
                      6 7 8 9 10 -1 -1 -1 -1 -1 -1)
    all_to_vir     = (-1 -1 -1 -1 0 1 2 3 4 5
                      -1 -1 -1 6 7 8 9 10 11)
  ********************************************************/

  // Orbital subspaces to MOs
  docc_to_mo.resize(ndocc);
  actv_to_mo.resize(nactv);
  extr_to_mo.resize(nextr);
  occ_to_mo.resize(nocc);
  vir_to_mo.resize(nvir);
  all_to_mo.resize(nall);

  int  index      = 0;
  int  mo_index   = 0;
  int  docc_index = 0;
  int  actv_index = 0;
  int  extr_index = 0;
  int  occ_index  = 0;
  int  vir_index  = 0;
  int  all_index  = 0;

  for(int h = 0; h < nirreps; ++h){   // Loop over irreps
    index += focc[h];                   // Offset by the focc

    for(int i = 0; i < docc[h]; ++i){   //Loop over doubly-occupied MOs
      docc_to_mo[docc_index] = index;
      occ_to_mo[occ_index]   = index;
      all_to_mo[all_index]   = index;
      docc_index++;
      occ_index++;
      all_index++;
      index++;
    }
    for(int i = 0; i < actv[h]; ++i){   //Loop over active MOs
      actv_to_mo[actv_index]=index;
      occ_to_mo[occ_index]  =index;
      vir_to_mo[vir_index]  =index;
      all_to_mo[all_index]  =index;
      actv_index++;
      occ_index++;
      vir_index++;
      all_index++;
      index++;
    }
    for(int i = 0; i < extr[h]; ++i){   //Loop over external MOs
      vir_to_mo[vir_index]  =index;
      extr_to_mo[extr_index]=index;
      all_to_mo[all_index]  =index;
      vir_index++;
      extr_index++;
      all_index++;
      index++;
    }
    index += fvir[h];                   // Offset by the frozen virtual MOs
  }



  // Orbital subspaces within themselves

  // Define the occ -> vir mapping
  occ_to_vir.resize(nocc,-1);
  for(int i = 0; i < nocc; ++i){
    for(int a = 0; a < nvir; ++a)
      if(occ_to_mo[i] == vir_to_mo[a])
        occ_to_vir[i] = a;
  }

  // Define the all -> occ and all -> vir mappings
  all_to_vir.resize(nall,-1);
  all_to_occ.resize(nall,-1);
  for(int p = 0; p < nall; ++p){
    for(int i = 0; i < nocc; ++i)
      if(all_to_mo[p]==occ_to_mo[i])
        all_to_occ[p] = i;
    for(int a = 0; a < nvir; ++a)
      if(all_to_mo[p]==vir_to_mo[a])
        all_to_vir[p] = a;
  }

  // Define the act -> occ and act -> vir mappings
  actv_to_occ.resize(nactv);
  actv_to_vir.resize(nactv);
  for(int u = 0; u < nactv; ++u){
    for(int i = 0; i < nocc; ++i)
      if(actv_to_mo[u]==occ_to_mo[i])
        actv_to_occ[u] = i;
    for(int a = 0; a < nvir; ++a)
      if(actv_to_mo[u]==vir_to_mo[a])
        actv_to_vir[u] = a;
  }

  // Define the occ -> act and vir -> act mappings
  occ_to_actv.resize(nocc,-1);
  vir_to_actv.resize(nvir,-1);
  for(int i = 0; i < nocc; ++i){
    for(int u = 0; u < nactv; ++u)
      if(occ_to_mo[i]==actv_to_mo[u])
        occ_to_actv[i]=u;
  }
  for(int a = 0; a < nvir; ++a){
    for(int u = 0; u < nactv; ++u)
      if(vir_to_mo[a]==actv_to_mo[u])
        vir_to_actv[a]=u;
  }

  // These will tell you if a certain orbital in the
  // generalized occupied space is active or not.
  // Used by the perturbation theory code.
  is_actv_in_occ.resize(nocc,false);
  is_actv_in_vir.resize(nvir,false);
  for(int u = 0; u < nactv; ++u){
    is_actv_in_occ[actv_to_occ[u]] = true;
    is_actv_in_vir[actv_to_vir[u]] = true;
  }
}

}

// OBSOLETE CODE
//first_occupied_mo[alpha] = new int[nirreps];
//first_active_mo[alpha]   = new int[nirreps];
//first_virtual_mo[alpha]  = new int[nirreps];
//last_occupied_mo[alpha]  = new int[nirreps];
//last_active_mo[alpha]    = new int[nirreps];
//last_virtual_mo[alpha]   = new int[nirreps];


/********************************************************
  Build the first and last arrays. These quantities will
  give the first and the last element of a mo space
  (doubly occupied, active, external, occupied, virtual)
  for a given irrep. For the example above:

  first_occupied_mo[alpha] = (4 18)
  first_active_mo[alpha]   = (8 21)
  first_virtual_mo[alpha]  = (8 21)
  last_occupied_mo[alpha]  = (10 23)
  last_active_mo[alpha]    = (10 23)
  last_virtual_mo[alpha]   = (14 27)
********************************************************/

//first_orbs_mo.resize(nirreps);
//first_so_mo.resize(nirreps);
//last_orbs_mo.resize(nirreps);
//last_so_mo.resize(nirreps);
//
//
//
//int sum = 0;
//for(int i=0;i<nirreps;i++){        //Loop over irreps
//  sum+=focc[i];
//  first_occupied_mo[alpha][i]=sum;
//  last_occupied_mo[alpha][i]=sum+docc[i]+actv[i];
//  sum+=docc[i];
//  first_active_mo[alpha][i]=sum;
//  last_active_mo[alpha][i]=sum+actv[i];
//  first_virtual_mo[alpha][i]=sum;
//  last_virtual_mo[alpha][i]= sum + actv[i] + extr[i];
//  sum += actv[i] + extr[i] + fvir[i];
//}
//sum = 0;
//for(int i=0;i<nirreps;i++){        //Loop over irreps
//  first_orbs_mo[i]=sum;
//  sum+=orbspi[i];
//  last_orbs_mo[i]=sum;
//}
//sum = 0;
//for(int i=0;i<nirreps;i++){        //Loop over irreps
//  first_so_mo[i]=sum;
//  sum+=sopi[i];
//  last_so_mo[i]=sum;
//}
//
//// Copy the alpha array to beta. This works only for MRCC
//for(int i=0;i<nirreps;i++){        //Loop over irreps
//  first_occupied_mo[beta] = first_occupied_mo[alpha];
//  first_active_mo[beta]   = first_active_mo[alpha];
//  first_virtual_mo[beta]  = first_virtual_mo[alpha];
//  last_occupied_mo[beta]  = last_occupied_mo[alpha];
//  last_active_mo[beta]    = last_active_mo[alpha];
//  last_virtual_mo[beta]   = last_virtual_mo[alpha];
//}

