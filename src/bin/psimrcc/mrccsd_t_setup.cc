#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>

#include "blas.h"
#include "heff.h"
#include "index_iterator.h"
#include "matrix.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

namespace psi{ namespace psimrcc{

void MRCCSD_T::startup()
{
  if(options_get_str("TRIPLES_ALGORITHM") == "SPIN_ADAPTED"){
    triples_algorithm = SpinAdaptedTriples;
  }else if(options_get_str("TRIPLES_ALGORITHM") == "RESTRICTED"){
    triples_algorithm = RestrictedTriples;
  }else{
    triples_algorithm = UnrestrictedTriples;
  }

  nirreps   = moinfo->get_nirreps();
  nrefs     = moinfo->get_ref_size(AllRefs);
  threshold = 0.1 * pow(10.0,-static_cast<double>(options_get_int("CONVERGENCE")));

  build_W_intermediates();

  o   = blas->get_index("[o]");
  oo  = blas->get_index("[oo]");
  ov  = blas->get_index("[ov]");
  v   = blas->get_index("[v]");
  vo  = blas->get_index("[vo]");
  vv  = blas->get_index("[vv]");
  vvv = blas->get_index("[vvv]");
  ovv = blas->get_index("[ovv]");
  ooo = blas->get_index("[ooo]");

  T2_ij_a_b = new IndexMatrix();
  T2_iJ_a_B = new IndexMatrix();
  T2_iJ_B_a = new IndexMatrix();
  T2_IJ_A_B = new IndexMatrix();

  T2_i_ab_j = new IndexMatrix();
  T2_i_aB_J = new IndexMatrix();
  T2_J_aB_i = new IndexMatrix();
  T2_I_AB_J = new IndexMatrix();

  V_k_bc_e = new IndexMatrix();
  V_k_bC_E = new IndexMatrix();
  V_K_bC_e = new IndexMatrix();

  V_jk_c_m = new IndexMatrix();
  V_jK_c_M = new IndexMatrix();
  V_jK_C_m = new IndexMatrix();

  form_T2_ij_a_b(T2_ij_a_b,true,true,false);
  form_T2_ij_a_b(T2_iJ_a_B,true,false,false);
  form_T2_ij_a_b(T2_iJ_B_a,true,false,true);   // T2_iJ_B_a = t2_{iJ}^{aB}
  form_T2_ij_a_b(T2_IJ_A_B,false,false,false);

  form_T2_i_ab_j(T2_i_ab_j,true,true,false);   // T2_i_ab_j = t2_{ij}^{ab}
  form_T2_i_ab_j(T2_i_aB_J,true,false,false);  // T2_i_aB_J = t2_{iJ}^{aB}
  form_T2_i_ab_j(T2_J_aB_i,true,false,true);   // T2_J_aB_i = t2_{iJ}^{aB}
  form_T2_i_ab_j(T2_I_AB_J,false,false,false); // T2_I_AB_J = t2_{IJ}^{AB}

  form_V_k_bc_e(V_k_bc_e,1.0,-1.0); // = <bc:ek>
  form_V_k_bc_e(V_k_bC_E,0.0,1.0);  // = <bC|kE> = <Ek|Cb>
  form_V_k_bc_e(V_K_bC_e,1.0,0.0);  // = <bC|eK> = <eK|bC>

  form_V_jk_c_m(V_jk_c_m,1.0,-1.0); // = <jk:mc>
  form_V_jk_c_m(V_jK_c_M,0.0,1.0);  // = <jk|cm>
  form_V_jk_c_m(V_jK_C_m,1.0,0.0);  // = <jk|mc>

  for(int mu = 0; mu < nrefs; ++mu){
    int unique_mu = moinfo->get_ref_number(mu,AllRefs);

    //  Unique references
    if(mu == unique_mu){
      // Setup the denominators
      double*** F_oo = blas->get_MatTmp("fock[o][o]",mu,none)->get_matrix();
      std::vector<double>  e_oo_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_oo_mu.push_back(F_oo[i_sym][i_rel][i_rel]);
        }
      }
      e_oo.push_back(e_oo_mu);

      double*** F_OO = blas->get_MatTmp("fock[O][O]",mu,none)->get_matrix();
      std::vector<double>  e_OO_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_OO_mu.push_back(F_OO[i_sym][i_rel][i_rel]);
        }
      }
      e_OO.push_back(e_OO_mu);


      double*** F_vv = blas->get_MatTmp("fock[v][v]",mu,none)->get_matrix();
      std::vector<double>  e_vv_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_vv_mu.push_back(F_vv[a_sym][a_rel][a_rel]);
        }
      }
      e_vv.push_back(e_vv_mu);

      double*** F_VV = blas->get_MatTmp("fock[V][V]",mu,none)->get_matrix();
      std::vector<double>  e_VV_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_VV_mu.push_back(F_VV[a_sym][a_rel][a_rel]);
        }
      }
      e_VV.push_back(e_VV_mu);

      F_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
      F_OV.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());

      if(options_get_bool("HEFF4")){
        F2_ov.push_back(blas->get_MatTmp("F_me[o][v]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("F_ME[O][V]",unique_mu,none)->get_matrix());
      }else{
        F2_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
      }

      W_ooov.push_back(blas->get_MatTmp("W_ijka[oo][ov]",unique_mu,none)->get_matrix());
      W_oOoV.push_back(blas->get_MatTmp("W_iJkA[oO][oV]",unique_mu,none)->get_matrix());
      W_OoOv.push_back(blas->get_MatTmp("W_IjKa[Oo][Ov]",unique_mu,none)->get_matrix());
      W_OOOV.push_back(blas->get_MatTmp("W_IJKA[OO][OV]",unique_mu,none)->get_matrix());

      W_vovv.push_back(blas->get_MatTmp("W_aibc[v][ovv]",unique_mu,none)->get_matrix());
      W_vOvV.push_back(blas->get_MatTmp("W_aIbC[v][OvV]",unique_mu,none)->get_matrix());
      W_VoVv.push_back(blas->get_MatTmp("W_AiBc[V][oVv]",unique_mu,none)->get_matrix());
      W_VOVV.push_back(blas->get_MatTmp("W_AIBC[V][OVV]",unique_mu,none)->get_matrix());

      T1_ov.push_back(blas->get_MatTmp("t1[o][v]",mu,none)->get_matrix());
      T1_OV.push_back(blas->get_MatTmp("t1[O][V]",mu,none)->get_matrix());

      T2_oovv.push_back(blas->get_MatTmp("t2[oo][vv]",mu,none)->get_matrix());
      T2_oOvV.push_back(blas->get_MatTmp("t2[oO][vV]",mu,none)->get_matrix());
      T2_OOVV.push_back(blas->get_MatTmp("t2[OO][VV]",mu,none)->get_matrix());

      is_aocc.push_back(moinfo->get_is_aocc(mu,AllRefs));
      is_bocc.push_back(moinfo->get_is_bocc(mu,AllRefs));
      is_avir.push_back(moinfo->get_is_avir(mu,AllRefs));
      is_bvir.push_back(moinfo->get_is_bvir(mu,AllRefs));
    }else{
      // Setup the denominators
      double*** F_oo = blas->get_MatTmp("fock[O][O]",unique_mu,none)->get_matrix();
      std::vector<double>  e_oo_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_oo_mu.push_back(F_oo[i_sym][i_rel][i_rel]);
        }
      }
      e_oo.push_back(e_oo_mu);

      double*** F_OO = blas->get_MatTmp("fock[o][o]",unique_mu,none)->get_matrix();
      std::vector<double>  e_OO_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_OO_mu.push_back(F_OO[i_sym][i_rel][i_rel]);
        }
      }
      e_OO.push_back(e_OO_mu);

      double*** F_vv = blas->get_MatTmp("fock[V][V]",unique_mu,none)->get_matrix();
      std::vector<double>  e_vv_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_vv_mu.push_back(F_vv[a_sym][a_rel][a_rel]);
        }
      }
      e_vv.push_back(e_vv_mu);

      double*** F_VV = blas->get_MatTmp("fock[v][v]",unique_mu,none)->get_matrix();
      std::vector<double>  e_VV_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_VV_mu.push_back(F_VV[a_sym][a_rel][a_rel]);
        }
      }
      e_VV.push_back(e_VV_mu);

      F_ov.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
      F_OV.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());

      if(options_get_bool("HEFF4")){
        F2_ov.push_back(blas->get_MatTmp("F_ME[O][V]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("F_me[o][v]",unique_mu,none)->get_matrix());
      }else{
        F2_ov.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
      }

      W_ooov.push_back(blas->get_MatTmp("W_IJKA[OO][OV]",unique_mu,none)->get_matrix());
      W_oOoV.push_back(blas->get_MatTmp("W_IjKa[Oo][Ov]",unique_mu,none)->get_matrix());
      W_OoOv.push_back(blas->get_MatTmp("W_iJkA[oO][oV]",unique_mu,none)->get_matrix());
      W_OOOV.push_back(blas->get_MatTmp("W_ijka[oo][ov]",unique_mu,none)->get_matrix());

      W_vovv.push_back(blas->get_MatTmp("W_AIBC[V][OVV]",unique_mu,none)->get_matrix());
      W_vOvV.push_back(blas->get_MatTmp("W_AiBc[V][oVv]",unique_mu,none)->get_matrix());
      W_VoVv.push_back(blas->get_MatTmp("W_aIbC[v][OvV]",unique_mu,none)->get_matrix());
      W_VOVV.push_back(blas->get_MatTmp("W_aibc[v][ovv]",unique_mu,none)->get_matrix());

      T1_ov.push_back(blas->get_MatTmp("t1[O][V]",unique_mu,none)->get_matrix());
      T1_OV.push_back(blas->get_MatTmp("t1[o][v]",unique_mu,none)->get_matrix());

      T2_oovv.push_back(blas->get_MatTmp("t2[OO][VV]",unique_mu,none)->get_matrix());
      T2_oOvV.push_back(blas->get_MatTmp("t2[Oo][Vv]",unique_mu,none)->get_matrix());
      T2_OOVV.push_back(blas->get_MatTmp("t2[oo][vv]",unique_mu,none)->get_matrix());

      is_bocc.push_back(moinfo->get_is_aocc(unique_mu,AllRefs));
      is_aocc.push_back(moinfo->get_is_bocc(unique_mu,AllRefs));
      is_bvir.push_back(moinfo->get_is_avir(unique_mu,AllRefs));
      is_avir.push_back(moinfo->get_is_bvir(unique_mu,AllRefs));
    }

    if(options_get_str("CORR_CCSD_T") == "STANDARD"){
      vector<double> factor_row;
      for(int nu = 0; nu < nrefs; ++nu){
        double c_mu   = h_eff->get_right_eigenvector(mu);
        double c_nu   = h_eff->get_right_eigenvector(nu);
        double factor = h_eff->get_matrix(mu,nu) * c_nu / c_mu;
        if(options_get_bool("TIKHONOW_TRIPLES")){
          double omega  = static_cast<double>(options_get_int("TIKHONOW_OMEGA")) / 1000.0;
          factor = h_eff->get_matrix(mu,nu) * c_nu * c_mu / (pow(c_mu,2.0) + pow(omega,2.0));
        }
        factor_row.push_back(factor);
      }
      Mk_factor.push_back(factor_row);

      Mk_shift.push_back(h_eff->get_eigenvalue() - h_eff->get_matrix(mu,mu));
    }else if(options_get_str("CORR_CCSD_T") == "PITTNER"){
      vector<double> factor_row;
      for(int nu = 0; nu < nrefs; ++nu){
        factor_row.push_back(0.0);
      }
      Mk_factor.push_back(factor_row);

      Mk_shift.push_back(0.0);
    }

    vector<double> d_h_eff_row;
    for(int nu = 0; nu < nrefs; ++nu){
      d_h_eff_row.push_back(0.0);
    }
    d_h_eff.push_back(d_h_eff_row);
  }

  V_oovv = blas->get_MatTmp("<[oo]:[vv]>",none)->get_matrix();
  V_oOvV = blas->get_MatTmp("<[oo]|[vv]>",none)->get_matrix();

  // Allocate Z, this will hold the results
  allocate2(BlockMatrix**,Z,nrefs,nirreps);
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      Z[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate W
  if((triples_algorithm == UnrestrictedTriples) or (triples_algorithm == RestrictedTriples)){
    allocate2(BlockMatrix**,W,nrefs,nirreps);
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        W[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
      }
    }
  }else if(triples_algorithm == SpinAdaptedTriples){
    allocate2(BlockMatrix**,W_ijk,nrefs,nirreps);
    allocate2(BlockMatrix**,W_ikj,nrefs,nirreps);
    allocate2(BlockMatrix**,W_jki,nrefs,nirreps);
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        W_ijk[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
        W_ikj[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
        W_jki[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
      }
    }
  }

  // Allocate T
  allocate2(BlockMatrix**,T,nrefs,nirreps);
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      T[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  e4T.assign(nrefs,0.0);
  e4ST.assign(nrefs,0.0);
  e4DT.assign(nrefs,0.0);

  E4T_ooo.assign(nrefs,0.0);
  E4T_ooO.assign(nrefs,0.0);
  E4T_oOO.assign(nrefs,0.0);
  E4T_OOO.assign(nrefs,0.0);
  E4ST_ooo.assign(nrefs,0.0);
  E4ST_ooO.assign(nrefs,0.0);
  E4ST_oOO.assign(nrefs,0.0);
  E4ST_OOO.assign(nrefs,0.0);
  E4DT_ooo.assign(nrefs,0.0);
  E4DT_ooO.assign(nrefs,0.0);
  E4DT_oOO.assign(nrefs,0.0);
  E4DT_OOO.assign(nrefs,0.0);
  E4_ooo.assign(nrefs,0.0);
  E4_ooO.assign(nrefs,0.0);
  E4_oOO.assign(nrefs,0.0);
  E4_OOO.assign(nrefs,0.0);
}



void MRCCSD_T::build_W_intermediates()
{
  blas->solve("W_ijka[oo][ov]{u}  = <[oo]:[ov]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_ijka[oo][ov]{u} += #4123# <[v]:[voo]> 1@2 t1[o][v]{u}");

  blas->solve("W_iJkA[oO][oV]{u}  = <[oo]|[ov]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_iJkA[oO][oV]{u} += #4123# <[v]|[voo]> 1@2 t1[o][v]{u}");

  blas->solve("W_IjKa[Oo][Ov]{u}  = <[oo]|[ov]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_IjKa[Oo][Ov]{u} += #4123# <[v]|[voo]> 1@2 t1[O][V]{u}");

  blas->solve("W_IJKA[OO][OV]{u}  = <[oo]:[ov]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_IJKA[OO][OV]{u} += #4123# <[v]:[voo]> 1@2 t1[O][V]{u}");

  blas->solve("W_aibc[v][ovv]{u}  = <[v]:[ovv]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_aibc[v][ovv]{u} += - t1[o][v]{u} 1@1 <[o]:[ovv]>");

  blas->solve("W_aIbC[v][OvV]{u}  = <[v]|[ovv]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_aIbC[v][OvV]{u} += - t1[o][v]{u} 1@1 <[o]|[ovv]>");

  blas->solve("W_AiBc[V][oVv]{u}  = <[v]|[ovv]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_AiBc[V][oVv]{u} += - t1[O][V]{u} 1@1 <[o]|[ovv]>");

  blas->solve("W_AIBC[V][OVV]{u}  = <[v]:[ovv]>");
  if(options_get_bool("HEFF4"))
    blas->solve("W_AIBC[V][OVV]{u} += - t1[O][V]{u} 1@1 <[o]:[ovv]>");
}

void MRCCSD_T::cleanup()
{
  delete T2_ij_a_b;
  delete T2_iJ_a_B;
  delete T2_iJ_B_a;
  delete T2_IJ_A_B;

  delete T2_i_ab_j;
  delete T2_i_aB_J;
  delete T2_J_aB_i;
  delete T2_I_AB_J;

  delete V_k_bc_e;
  delete V_k_bC_E;
  delete V_K_bC_e;

  delete V_jk_c_m;
  delete V_jK_c_M;
  delete V_jK_C_m;

  // Deallocate Z
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete Z[mu][h];
    }
  }
  release2(Z);

  // Deallocate W
  if((triples_algorithm == UnrestrictedTriples) or (triples_algorithm == RestrictedTriples)){
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        delete W[mu][h];
      }
    }
    release2(W);
  }else if(triples_algorithm == SpinAdaptedTriples){
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        delete W_ijk[mu][h];  delete W_ikj[mu][h];  delete W_jki[mu][h];
      }
    }
    release2(W_ijk);  release2(W_ikj);  release2(W_jki);
  }

  // Deallocate T
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete T[mu][h];
    }
  }
  release2(T);
}

}} /* End Namespaces */
