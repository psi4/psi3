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
//  | focc  | docc  |act|    ext    |fvir |focc | docc   | act |    ext    | fvir|
//  |-------|     occ   |-----------------------|       occ    |-----------------|
//  |---------------|     vir       |--------------------|     vir         |-----|
//  |-------|          all          |-----------|            all           |-----|
//  |----------------------------------------------------------------------------|
//

  /********************************************************
    Build the first and last arrays. These quantities will
    give the first and the last element of a mo space
    (doubly occupied, active, external, occupied, virtual)
    for a given irrep. For the example above:

    first_occupied_pitzer[alpha] = (4 18)
    first_active_pitzer[alpha]   = (8 21)
    first_virtual_pitzer[alpha]  = (8 21)
    last_occupied_pitzer[alpha]  = (10 23)
    last_active_pitzer[alpha]    = (10 23)
    last_virtual_pitzer[alpha]   = (14 27)
  ********************************************************/

  first_orbs_pitzer            = new int[nirreps];
  first_so_pitzer              = new int[nirreps];
  first_occupied_pitzer[alpha] = new int[nirreps];
  first_active_pitzer[alpha]   = new int[nirreps];
  first_virtual_pitzer[alpha]  = new int[nirreps];
  last_orbs_pitzer             = new int[nirreps];
  last_so_pitzer               = new int[nirreps];
  last_occupied_pitzer[alpha]  = new int[nirreps];
  last_active_pitzer[alpha]    = new int[nirreps];
  last_virtual_pitzer[alpha]   = new int[nirreps];

  int sum = 0;
  for(int i=0;i<nirreps;i++){        //Loop over irreps
    sum+=focc[i];
    first_occupied_pitzer[alpha][i]=sum;
    last_occupied_pitzer[alpha][i]=sum+docc[i]+actv[i];
    sum+=docc[i];
    first_active_pitzer[alpha][i]=sum;
    last_active_pitzer[alpha][i]=sum+actv[i];
    first_virtual_pitzer[alpha][i]=sum;
    last_virtual_pitzer[alpha][i]=sum+actv[i]+avir[i];
    sum+=actv[i]+avir[i]+fvir[i];
  }
  sum = 0;
  for(int i=0;i<nirreps;i++){        //Loop over irreps
    first_orbs_pitzer[i]=sum;
    sum+=orbspi[i];
    last_orbs_pitzer[i]=sum;
  }
  sum = 0;
  for(int i=0;i<nirreps;i++){        //Loop over irreps
    first_so_pitzer[i]=sum;
    sum+=sopi[i];
    last_so_pitzer[i]=sum;
  }

  // Copy the alpha array to beta. This works only for MRCC
  for(int i=0;i<nirreps;i++){        //Loop over irreps
    first_occupied_pitzer[beta] = first_occupied_pitzer[alpha];
    first_active_pitzer[beta]   = first_active_pitzer[alpha];
    first_virtual_pitzer[beta]  = first_virtual_pitzer[alpha];
    last_occupied_pitzer[beta]  = last_occupied_pitzer[alpha];
    last_active_pitzer[beta]    = last_active_pitzer[alpha];
    last_virtual_pitzer[beta]   = last_virtual_pitzer[alpha];
  }

  /********************************************************
    Build the arrays that connect the subspaces of the
    non-frozen MOs to the the non-frozen list of MOs
    in Pitzer. This is what you need in all computations.

    For the example above:

    docc_to_pitzer = (4 5 6 7 18 19 20)
    act_to_pitzer  = (8 9 21 22)
    ext_to_pitzer  = (10 11 12 13 23 24 25 26)
    occ_to_pitzer  = (4 5 6 7 8 9 18 19 20 21 22)
    vir_to_pitzer  = (8 9 10 11 12 13 21 22 23 24 25 26)
    all_to_pitzer  = (4 5 6 7 8 9 10 11 12 13 18 19 20
                      21 22 23 24 25 26)
    all_to_occ     = (0 1 2 3 4 5 -1 -1 -1 -1
                      6 7 8 9 10 -1 -1 -1 -1 -1 -1)
    all_to_vir     = (-1 -1 -1 -1 0 1 2 3 4 5
                      -1 -1 -1 6 7 8 9 10 11)
  ********************************************************/

  so_to_pitzer   = new int[nso];
  orbs_to_pitzer = new int[norbs];
  docc_to_pitzer = new int[ndocc];
  act_to_pitzer = new int[nactv];
  ext_to_pitzer = new int[navir];
  occ_to_pitzer = new int[nocc];
  vir_to_pitzer = new int[nvir];
  all_to_pitzer = new int[nmo];
  all_to_vir    = new int[nmo];
  all_to_occ    = new int[nmo];
  occ_to_vir    = new int[nocc];

  act_to_occ    = new int[nactv];
  act_to_vir    = new int[nactv];
  occ_to_act    = new int[nocc];
  vir_to_act    = new int[nvir];

  for(int i=0;i<nso;i++)
    so_to_pitzer[i]=i;

  sum=0;
  int  orbs_index = 0;
  int  docc_index = 0;
  int  act_index = 0;
  int  ext_index = 0;
  int  occ_index = 0;
  int  vir_index = 0;
  int  all_index = 0;

  for(int i=0;i<nirreps;i++){        //Loop over irreps
    for(int j=0;j<focc[i];j++){              //Loop over focc
      orbs_to_pitzer[orbs_index]=sum;
      orbs_index++;
      sum++;
    }
    for(int j=0;j<docc[i];j++){              //Loop over docc
      docc_to_pitzer[docc_index]=sum;
      occ_to_pitzer[occ_index]=sum;
      all_to_pitzer[all_index]=sum;
      orbs_to_pitzer[orbs_index]=sum;
      orbs_index++;
      docc_index++;
      occ_index++;
      all_index++;
      sum++;
    }
    for(int j=0;j<actv[i];j++){              //Loop over act
      act_to_pitzer[act_index]=sum;
      occ_to_pitzer[occ_index]=sum;
      vir_to_pitzer[vir_index]=sum;
      all_to_pitzer[all_index]=sum;
      orbs_to_pitzer[orbs_index]=sum;
      orbs_index++;
      act_index++;
      occ_index++;
      vir_index++;
      all_index++;
      sum++;
    }
    for(int j=0;j<avir[i];j++){              //Loop over ext
      vir_to_pitzer[vir_index]=sum;
      ext_to_pitzer[ext_index]=sum;
      all_to_pitzer[all_index]=sum;
      orbs_to_pitzer[orbs_index]=sum;
      orbs_index++;
      vir_index++;
      ext_index++;
      all_index++;
      sum++;
    }
    for(int j=0;j<fvir[i];j++){              //Loop over fvir
      orbs_to_pitzer[orbs_index]=sum;
      orbs_index++;
      sum++;
    }
  }

  for(int i=0;i<nocc;i++){
    occ_to_vir[i]=-1;
    for(int j=0;j<nvir;j++)
      if(occ_to_pitzer[i]==vir_to_pitzer[j])
        occ_to_vir[i]=j;
  }

  // Define the all_to_occ and all_to_vir mappings
  for(int i=0;i<nmo;i++){
    all_to_occ[i]=-1;
    all_to_vir[i]=-1;
    for(int j=0;j<nocc;j++)
      if(all_to_pitzer[i]==occ_to_pitzer[j])
        all_to_occ[i]=j;
    for(int j=0;j<nvir;j++)
      if(all_to_pitzer[i]==vir_to_pitzer[j])
        all_to_vir[i]=j;
  }
  
  // Define the act_to_occ and act_to_vir mappings
  for(int i=0;i<nactv;i++){
    for(int j=0;j<nocc;j++)
      if(act_to_pitzer[i]==occ_to_pitzer[j])
        act_to_occ[i]=j;
    for(int j=0;j<nvir;j++)
      if(act_to_pitzer[i]==vir_to_pitzer[j])
        act_to_vir[i]=j;
  }

  // Define the occ_to_act and vir_to_act mappings
  for(int i=0;i<nocc;i++){
    occ_to_act[i] = -1;
    for(int j=0;j<nactv;j++)
      if(occ_to_pitzer[i]==act_to_pitzer[j])
        occ_to_act[i]=j;
  }
  for(int i=0;i<nvir;i++){
    vir_to_act[i] = -1;
    for(int j=0;j<nactv;j++)
      if(vir_to_pitzer[i]==act_to_pitzer[j])
        vir_to_act[i]=j;
  }

  // These will tell you if a certain orbital in the
  // generalized occupied space is active or not.
  // Used by the perturbation theory code.
  is_act_in_occ = new bool[nocc];
  is_act_in_vir = new bool[nvir];
  for(int i=0;i<nocc;i++)
    is_act_in_occ[i]=false;
  for(int i=0;i<nvir;i++)
    is_act_in_vir[i]=false;
  for(int i=0;i<nactv;i++){
    is_act_in_occ[act_to_occ[i]]=true;
    is_act_in_vir[act_to_vir[i]]=true;
  }
}

}
