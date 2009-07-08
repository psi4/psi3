/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <cstdlib>

#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "debugging.h"
#include "index_iterator.h"
#include "mrcc.h"
#include "special_matrices.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

void CCMRCC::compute_perturbative_triples()
{
  Timer timer;
  fprintf(outfile,"\n\n  Computing (T) correction");

  //  Algorithm
  compute_ooo_triples();
  compute_OOO_triples();
  compute_ooO_triples();

  fprintf(outfile,"\n  Timing for triples: %20.6f s",timer.get());
  fflush(outfile);
}

void CCMRCC::compute_ooo_triples()
{
  IndexMatrix T2_ij_a_b;
  IndexMatrix T2_i_ab_j;
  IndexMatrix V_k_bc_e;
  IndexMatrix V_jk_c_m;

  form_T2_ij_a_b(&T2_ij_a_b,true,true,false);
  form_T2_i_ab_j(&T2_i_ab_j,true,true,false);
  form_V_k_bc_e(&V_k_bc_e,1.0,-1.0);
  form_V_jk_c_m(&V_jk_c_m,1.0,-1.0);

  int nirreps = moinfo->get_nirreps();
  int nurefs  = moinfo->get_nunique();

  CCIndex* o   = blas->get_index("[o]");
  CCIndex* oo  = blas->get_index("[oo]");
  CCIndex* v   = blas->get_index("[v]");
  CCIndex* vv  = blas->get_index("[vv]");
  CCIndex* vvv = blas->get_index("[vvv]");


  CCIndexIterator  ijk("[ooo]");
  CCIndexIterator  abc("[vvv]");

  // Allocate Z, this will hold the results
  BlockMatrix*** Z;
  allocate2(BlockMatrix**,Z,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      Z[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate W
  BlockMatrix*** W;
  allocate2(BlockMatrix**,W,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      W[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate T
  BlockMatrix*** T;
  allocate2(BlockMatrix**,T,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      T[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

//  std::vector<bool> is_aocc = moinfo->get_is_aocc(reference,UniqueRefs);
//  std::vector<bool> is_avir = moinfo->get_is_avir(reference,UniqueRefs);

  // TODO: Work on the denominator (also generalized denominators)
  // TODO: Use blas for matrix multiplications (perhaps call it contract)

  vector<vector<bool> > is_aocc;
  vector<vector<bool> > is_avir;

  vector<double***> F_oo;
  vector<double***> F_ov;
  vector<double***> F_vv;
  vector<double***> T1_ov;
  vector<double***> V_oovv;
  vector<double***> T2_oovv;


  for(int mu = 0; mu < nurefs; ++mu){
    int unique_mu = moinfo->get_ref_number(mu,UniqueRefs);

    V_oovv.push_back(blas->get_MatTmp("<[oo]:[vv]>",none)->get_matrix());

    F_oo.push_back(blas->get_MatTmp("fock[o][o]",unique_mu,none)->get_matrix());
    F_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
    F_vv.push_back(blas->get_MatTmp("fock[v][v]",unique_mu,none)->get_matrix());
    T1_ov.push_back(blas->get_MatTmp("t1[o][v]",unique_mu,none)->get_matrix());
    T2_oovv.push_back(blas->get_MatTmp("t2[oo][vv]",unique_mu,none)->get_matrix());

    is_aocc.push_back(moinfo->get_is_aocc(mu,UniqueRefs));
    is_avir.push_back(moinfo->get_is_avir(mu,UniqueRefs));
  }

  vector<double> E4T(nurefs,0.0);
  vector<double> E4ST(nurefs,0.0);
  vector<double> E4DT(nurefs,0.0);

  while(++ijk){
    int i_sym     = o->get_tuple_irrep(ijk.ind_abs[0]);
    int j_sym     = o->get_tuple_irrep(ijk.ind_abs[1]);
    int k_sym     = o->get_tuple_irrep(ijk.ind_abs[2]);

    size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs[0]);
    size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs[1]);
    size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs[2]);

    size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs[0]);
    size_t j_rel = o->get_tuple_rel_index(ijk.ind_abs[1]);
    size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs[2]);

    size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[1]);
    size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs[2],ijk.ind_abs[1]);
    size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[2]);

    int    jk_sym = oo->get_tuple_irrep(ijk.ind_abs[1],ijk.ind_abs[2]);
    size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs[1],ijk.ind_abs[2]);

    for(int mu = 0; mu < nurefs; ++mu){

      // Check if ijk belong to the occupied space of mu
      if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_aocc[mu][k_abs]){
        double D_ijk = F_oo[mu][i_sym][i_rel][i_rel]
                     + F_oo[mu][j_sym][j_rel][j_rel]
                     + F_oo[mu][k_sym][k_rel][k_rel];

        vector<bool>& is_aocc_mu = is_aocc[mu];
        vector<bool>& is_avir_mu = is_avir[mu];
        double***  F_ov_mu    = F_ov[mu];
        double***  F_vv_mu    = F_vv[mu];
        double***  T1_ov_mu   = T1_ov[mu];
        double***  V_oovv_mu  = V_oovv[mu];
        double***  T2_oovv_mu = T2_oovv[mu];

        // Compute P(k/ij) Z_ijk^abc
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(ij_abs,mu),
            V_k_bc_e.get_block_matrix(k_abs),0.0,1.0);
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(kj_abs,mu),
            V_k_bc_e.get_block_matrix(i_abs),1.0,-1.0);
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(ik_abs,mu),
            V_k_bc_e.get_block_matrix(j_abs),1.0,-1.0);

        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(ij_abs,mu),
            T2_i_ab_j.get_block_matrix(k_abs),1.0,-1.0);
        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(kj_abs,mu),
            T2_i_ab_j.get_block_matrix(i_abs),1.0,1.0);
        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(ik_abs,mu),
            T2_i_ab_j.get_block_matrix(j_abs),1.0,1.0);

        // Compute P(k/ij) Z_ijk^abc
        W[mu][ijk.sym]->cyclical_permutation_1_2(Z[mu][ijk.sym],vvv,v,vv);

        // Compute T
        abc.reset();
        abc.set_irrep(ijk.sym);
        while(++abc){
          int a_sym    = v->get_tuple_irrep(abc.ind_abs[0]);
          int b_sym    = v->get_tuple_irrep(abc.ind_abs[1]);
          int c_sym    = v->get_tuple_irrep(abc.ind_abs[2]);

          size_t a_abs = v->get_tuple_abs_index(abc.ind_abs[0]);
          size_t b_abs = v->get_tuple_abs_index(abc.ind_abs[1]);
          size_t c_abs = v->get_tuple_abs_index(abc.ind_abs[2]);

          size_t a_rel = v->get_tuple_rel_index(abc.ind_abs[0]);
          size_t b_rel = v->get_tuple_rel_index(abc.ind_abs[1]);
          size_t c_rel = v->get_tuple_rel_index(abc.ind_abs[2]);

          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);

          if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_avir_mu[c_abs]){
            double D_abc = F_vv_mu[a_sym][a_rel][a_rel]
                         + F_vv_mu[b_sym][b_rel][b_rel]
                         + F_vv_mu[c_sym][c_rel][c_rel];

            T[mu][ijk.sym]->set(a_sym,a_rel,bc_rel,W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel)/(D_ijk - D_abc));
          }
        }

        // Compute the Energy
        abc.reset();
        abc.set_irrep(ijk.sym);
        while(++abc){
          int a_sym     = v->get_tuple_irrep(abc.ind_abs[0]);

          size_t a_abs  = v->get_tuple_abs_index(abc.ind_abs[0]);
          size_t b_abs  = v->get_tuple_abs_index(abc.ind_abs[1]);
          size_t c_abs  = v->get_tuple_abs_index(abc.ind_abs[2]);

          size_t a_rel  = v->get_tuple_rel_index(abc.ind_abs[0]);

          int bc_sym = vv->get_tuple_irrep(abc.ind_abs[1],abc.ind_abs[2]);
          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);

          if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_avir_mu[c_abs]){
            // E4T
            E4T[mu]  += W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) *  T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) / 36.0;

            // E4ST
            if((i_sym == a_sym) & (jk_sym == bc_sym)){
              E4ST[mu] += 0.25 * T1_ov_mu[i_sym][i_rel][a_rel]
                        * V_oovv_mu[jk_sym][jk_rel][bc_rel]
                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
            }
            // E4DT
            if((i_sym == a_sym) & (jk_sym == bc_sym)){
              E4DT[mu] += 0.25 * F_ov_mu[i_sym][i_rel][a_rel]
                        * T2_oovv_mu[jk_sym][jk_rel][bc_rel]
                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
            }
          }
        }
      }
    }
  }

  double E4T_sum = 0.0;
  double E4ST_sum = 0.0;
  double E4DT_sum = 0.0;
  for(int mu = 0; mu < nurefs; ++mu){
    fprintf(outfile,"\n  E4T(%d)  = %20.15f",mu,E4T[mu]);
    fprintf(outfile,"\n  E4ST(%d) = %20.15f",mu,E4ST[mu]);
    fprintf(outfile,"\n  E4DT(%d) = %20.15f",mu,E4DT[mu]);
    E4T_sum += E4T[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4ST_sum += E4ST[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4DT_sum += E4DT[mu] * left_eigenvector[mu] * right_eigenvector[mu];
  }

  // Deallocate Z
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete Z[mu][h];
    }
  }
  release2(Z);

  // Deallocate W
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete W[mu][h];
    }
  }
  release2(W);

  // Deallocate T
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete T[mu][h];
    }
  }
  release2(T);
}

void CCMRCC::compute_OOO_triples()
{
  IndexMatrix T2_ij_a_b;
  IndexMatrix T2_i_ab_j;
  IndexMatrix V_k_bc_e;
  IndexMatrix V_jk_c_m;

  form_T2_ij_a_b(&T2_ij_a_b,false,false,false);
  form_T2_i_ab_j(&T2_i_ab_j,false,false,false);
  form_V_k_bc_e(&V_k_bc_e,1.0,-1.0);
  form_V_jk_c_m(&V_jk_c_m,1.0,-1.0);

  int nirreps = moinfo->get_nirreps();
  int nurefs  = moinfo->get_nunique();

  CCIndex* o   = blas->get_index("[o]");
  CCIndex* oo  = blas->get_index("[oo]");
  CCIndex* v   = blas->get_index("[v]");
  CCIndex* vv  = blas->get_index("[vv]");
  CCIndex* vvv = blas->get_index("[vvv]");


  CCIndexIterator  ijk("[ooo]");
  CCIndexIterator  abc("[vvv]");

  // Allocate Z, this will hold the results
  BlockMatrix*** Z;
  allocate2(BlockMatrix**,Z,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      Z[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate W
  BlockMatrix*** W;
  allocate2(BlockMatrix**,W,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      W[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate T
  BlockMatrix*** T;
  allocate2(BlockMatrix**,T,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      T[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

//  std::vector<bool> is_aocc = moinfo->get_is_aocc(reference,UniqueRefs);
//  std::vector<bool> is_avir = moinfo->get_is_avir(reference,UniqueRefs);

  // TODO: Work on the denominator (also generalized denominators)
  // TODO: Use blas for matrix multiplications (perhaps call it contract)

  vector<vector<bool> > is_bocc;
  vector<vector<bool> > is_bvir;

  vector<double***> F_OO;
  vector<double***> F_OV;
  vector<double***> F_VV;
  vector<double***> T1_OV;
  vector<double***> V_OOVV;
  vector<double***> T2_OOVV;


  for(int mu = 0; mu < nurefs; ++mu){
    int unique_mu = moinfo->get_ref_number(mu,UniqueRefs);

    V_OOVV.push_back(blas->get_MatTmp("<[oo]:[vv]>",none)->get_matrix());

    F_OO.push_back(blas->get_MatTmp("fock[O][O]",unique_mu,none)->get_matrix());
    F_OV.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
    F_VV.push_back(blas->get_MatTmp("fock[V][V]",unique_mu,none)->get_matrix());
    T1_OV.push_back(blas->get_MatTmp("t1[O][V]",unique_mu,none)->get_matrix());
    T2_OOVV.push_back(blas->get_MatTmp("t2[OO][VV]",unique_mu,none)->get_matrix());

    is_bocc.push_back(moinfo->get_is_bocc(mu,UniqueRefs));
    is_bvir.push_back(moinfo->get_is_bvir(mu,UniqueRefs));
  }

  vector<double> E4T(nurefs,0.0);
  vector<double> E4ST(nurefs,0.0);
  vector<double> E4DT(nurefs,0.0);

  while(++ijk){
    int i_sym     = o->get_tuple_irrep(ijk.ind_abs[0]);
    int j_sym     = o->get_tuple_irrep(ijk.ind_abs[1]);
    int k_sym     = o->get_tuple_irrep(ijk.ind_abs[2]);

    size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs[0]);
    size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs[1]);
    size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs[2]);

    size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs[0]);
    size_t j_rel = o->get_tuple_rel_index(ijk.ind_abs[1]);
    size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs[2]);

    size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[1]);
    size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs[2],ijk.ind_abs[1]);
    size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[2]);

    int    jk_sym = oo->get_tuple_irrep(ijk.ind_abs[1],ijk.ind_abs[2]);
    size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs[1],ijk.ind_abs[2]);

    for(int mu = 0; mu < nurefs; ++mu){

      // Check if ijk belong to the occupied space of mu
      if(is_bocc[mu][i_abs] && is_bocc[mu][j_abs] && is_bocc[mu][k_abs]){
        double D_ijk = F_OO[mu][i_sym][i_rel][i_rel]
                     + F_OO[mu][j_sym][j_rel][j_rel]
                     + F_OO[mu][k_sym][k_rel][k_rel];

        vector<bool>& is_bocc_mu = is_bocc[mu];
        vector<bool>& is_bvir_mu = is_bvir[mu];
        double***  F_OV_mu    = F_OV[mu];
        double***  F_VV_mu    = F_VV[mu];
        double***  T1_OV_mu   = T1_OV[mu];
        double***  V_OOVV_mu  = V_OOVV[mu];
        double***  T2_OOVV_mu = T2_OOVV[mu];

        // Compute P(k/ij) Z_ijk^abc
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(ij_abs,mu),
            V_k_bc_e.get_block_matrix(k_abs),0.0,1.0);
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(kj_abs,mu),
            V_k_bc_e.get_block_matrix(i_abs),1.0,-1.0);
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(ik_abs,mu),
            V_k_bc_e.get_block_matrix(j_abs),1.0,-1.0);

        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(ij_abs,mu),
            T2_i_ab_j.get_block_matrix(k_abs),1.0,-1.0);
        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(kj_abs,mu),
            T2_i_ab_j.get_block_matrix(i_abs),1.0,1.0);
        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(ik_abs,mu),
            T2_i_ab_j.get_block_matrix(j_abs),1.0,1.0);

        // Compute P(k/ij) Z_ijk^abc
        W[mu][ijk.sym]->cyclical_permutation_1_2(Z[mu][ijk.sym],vvv,v,vv);

        // Compute T
        abc.reset();
        abc.set_irrep(ijk.sym);
        while(++abc){
          int a_sym    = v->get_tuple_irrep(abc.ind_abs[0]);
          int b_sym    = v->get_tuple_irrep(abc.ind_abs[1]);
          int c_sym    = v->get_tuple_irrep(abc.ind_abs[2]);

          size_t a_abs = v->get_tuple_abs_index(abc.ind_abs[0]);
          size_t b_abs = v->get_tuple_abs_index(abc.ind_abs[1]);
          size_t c_abs = v->get_tuple_abs_index(abc.ind_abs[2]);

          size_t a_rel = v->get_tuple_rel_index(abc.ind_abs[0]);
          size_t b_rel = v->get_tuple_rel_index(abc.ind_abs[1]);
          size_t c_rel = v->get_tuple_rel_index(abc.ind_abs[2]);

          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);

          if(is_bvir_mu[a_abs] && is_bvir_mu[b_abs] && is_bvir_mu[c_abs]){
            double D_abc = F_VV_mu[a_sym][a_rel][a_rel]
                         + F_VV_mu[b_sym][b_rel][b_rel]
                         + F_VV_mu[c_sym][c_rel][c_rel];

            T[mu][ijk.sym]->set(a_sym,a_rel,bc_rel,W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel)/(D_ijk - D_abc));
          }
        }

        // Compute the Energy
        abc.reset();
        abc.set_irrep(ijk.sym);
        while(++abc){
          int a_sym     = v->get_tuple_irrep(abc.ind_abs[0]);

          size_t a_abs  = v->get_tuple_abs_index(abc.ind_abs[0]);
          size_t b_abs  = v->get_tuple_abs_index(abc.ind_abs[1]);
          size_t c_abs  = v->get_tuple_abs_index(abc.ind_abs[2]);

          size_t a_rel  = v->get_tuple_rel_index(abc.ind_abs[0]);

          int bc_sym = vv->get_tuple_irrep(abc.ind_abs[1],abc.ind_abs[2]);
          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);

          if(is_bvir_mu[a_abs] && is_bvir_mu[b_abs] && is_bvir_mu[c_abs]){
            // E4T
            E4T[mu]  += W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) *  T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) / 36.0;

            // E4ST
            if((i_sym == a_sym) & (jk_sym == bc_sym)){
              E4ST[mu] += 0.25 * T1_OV_mu[i_sym][i_rel][a_rel]
                        * V_OOVV_mu[jk_sym][jk_rel][bc_rel]
                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
            }
            // E4DT
            if((i_sym == a_sym) & (jk_sym == bc_sym)){
              E4DT[mu] += 0.25 * F_OV_mu[i_sym][i_rel][a_rel]
                        * T2_OOVV_mu[jk_sym][jk_rel][bc_rel]
                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
            }
          }
        }
      }
    }
  }

  double E4T_sum = 0.0;
  double E4ST_sum = 0.0;
  double E4DT_sum = 0.0;
  for(int mu = 0; mu < nurefs; ++mu){
    fprintf(outfile,"\n  E4T(%d)  = %20.15f",mu,E4T[mu]);
    fprintf(outfile,"\n  E4ST(%d) = %20.15f",mu,E4ST[mu]);
    fprintf(outfile,"\n  E4DT(%d) = %20.15f",mu,E4DT[mu]);
    E4T_sum += E4T[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4ST_sum += E4ST[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4DT_sum += E4DT[mu] * left_eigenvector[mu] * right_eigenvector[mu];
  }


  // Deallocate Z
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete Z[mu][h];
    }
  }
  release2(Z);

  // Deallocate W
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete W[mu][h];
    }
  }
  release2(W);

  // Deallocate T
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete T[mu][h];
    }
  }
  release2(T);
}


void CCMRCC::compute_ooO_triples()
{
  IndexMatrix T2_ij_a_b;
  IndexMatrix T2_iJ_a_B;
  IndexMatrix T2_iJ_B_a;

  IndexMatrix T2_J_aB_i;
  IndexMatrix T2_i_aB_J;


  IndexMatrix V_k_bc_e;
  IndexMatrix V_K_bC_e;
  IndexMatrix V_k_bC_E;

  IndexMatrix V_jk_c_m;
  IndexMatrix V_jK_c_M;

//  IndexMatrix V_jk_c_m;

  form_T2_ij_a_b(&T2_ij_a_b,true,true,false);
  form_T2_ij_a_b(&T2_iJ_a_B,true,false,false);
  form_T2_ij_a_b(&T2_iJ_B_a,true,false,true);  // T2_iJ_B_a = t2_{iJ}^{aB}

  form_T2_i_ab_j(&T2_J_aB_i,true,false,true);  // T2_J_aB_i = t2_{iJ}^{aB}
  form_T2_i_ab_j(&T2_i_aB_J,true,false,false);  // T2_i_aB_J = t2_{iJ}^{aB}

 // form_T2_i_ab_j(&T2_i_ab_j,true,true);

  form_V_k_bc_e(&V_k_bc_e,1.0,-1.0); // = <bc:ek>
  form_V_k_bc_e(&V_K_bC_e,1.0,0.0);  // = <bC|eK> = <eK|bC>
  form_V_k_bc_e(&V_k_bC_E,0.0,1.0);  // = <bC|kE> = <Ek|Cb>

  form_V_jk_c_m(&V_jk_c_m,1.0,-1.0);
  form_V_jk_c_m(&V_jK_c_M,0.0,1.0);
//
//  form_V_k_bc_e(&V_k_bc_e);
//  form_V_jk_c_m(&V_jk_c_m);

  int nirreps = moinfo->get_nirreps();
  int nurefs  = moinfo->get_nunique();

  CCIndex* o   = blas->get_index("[o]");
  CCIndex* oo  = blas->get_index("[oo]");
  CCIndex* v   = blas->get_index("[v]");
  CCIndex* vv  = blas->get_index("[vv]");
  CCIndex* vvv = blas->get_index("[vvv]");


  CCIndexIterator  ijk("[ooo]");
  CCIndexIterator  abc("[vvv]");

  // Allocate Z, this will hold the results
  BlockMatrix*** Z;
  allocate2(BlockMatrix**,Z,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      Z[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate W
  BlockMatrix*** W;
  allocate2(BlockMatrix**,W,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      W[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate T
  BlockMatrix*** T;
  allocate2(BlockMatrix**,T,nurefs,nirreps);
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      T[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

//  std::vector<bool> is_aocc = moinfo->get_is_aocc(reference,UniqueRefs);
//  std::vector<bool> is_avir = moinfo->get_is_avir(reference,UniqueRefs);

  // TODO: Work on the denominator (also generalized denominators)
  // TODO: Use blas for matrix multiplications (perhaps call it contract)

  vector<vector<bool> > is_aocc;
  vector<vector<bool> > is_avir;
  vector<vector<bool> > is_bocc;
  vector<vector<bool> > is_bvir;

  vector<double***> F_oo;
  vector<double***> F_ov;
  vector<double***> F_vv;
  vector<double***> T1_ov;
  vector<double***> V_oovv;
  vector<double***> T2_oovv;

  vector<double***> F_OO;
  vector<double***> F_OV;
  vector<double***> F_VV;
  vector<double***> T1_OV;
  vector<double***> V_OOVV;
  vector<double***> T2_OOVV;

  for(int mu = 0; mu < nurefs; ++mu){
    int unique_mu = moinfo->get_ref_number(mu,UniqueRefs);

    V_oovv.push_back(blas->get_MatTmp("<[oo]:[vv]>",none)->get_matrix());

    F_oo.push_back(blas->get_MatTmp("fock[o][o]",unique_mu,none)->get_matrix());
    F_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
    F_vv.push_back(blas->get_MatTmp("fock[v][v]",unique_mu,none)->get_matrix());
    T1_ov.push_back(blas->get_MatTmp("t1[o][v]",unique_mu,none)->get_matrix());
    T2_oovv.push_back(blas->get_MatTmp("t2[oo][vv]",unique_mu,none)->get_matrix());

    is_aocc.push_back(moinfo->get_is_aocc(mu,UniqueRefs));
    is_avir.push_back(moinfo->get_is_avir(mu,UniqueRefs));
    is_bocc.push_back(moinfo->get_is_bocc(mu,UniqueRefs));
    is_bvir.push_back(moinfo->get_is_bvir(mu,UniqueRefs));
  }

  vector<double> E4T(nurefs,0.0);
  vector<double> E4ST(nurefs,0.0);
  vector<double> E4DT(nurefs,0.0);

  while(++ijk){
    int i_sym     = o->get_tuple_irrep(ijk.ind_abs[0]);
    int j_sym     = o->get_tuple_irrep(ijk.ind_abs[1]);
    int k_sym     = o->get_tuple_irrep(ijk.ind_abs[2]);

    size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs[0]);
    size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs[1]);
    size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs[2]);

    size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs[0]);
    size_t j_rel = o->get_tuple_rel_index(ijk.ind_abs[1]);
    size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs[2]);

    size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[1]);
    size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs[2],ijk.ind_abs[1]);
    size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs[0],ijk.ind_abs[2]);
    size_t jk_abs = oo->get_tuple_abs_index(ijk.ind_abs[1],ijk.ind_abs[2]);

    int    jk_sym = oo->get_tuple_irrep(ijk.ind_abs[1],ijk.ind_abs[2]);
    size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs[1],ijk.ind_abs[2]);

    for(int mu = 0; mu < nurefs; ++mu){

      // Check if ijk belong to the occupied space of mu
      if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]){
        double D_ijK = F_oo[mu][i_sym][i_rel][i_rel]
                     + F_oo[mu][j_sym][j_rel][j_rel]
                     + F_oo[mu][k_sym][k_rel][k_rel];

        vector<bool>& is_aocc_mu = is_aocc[mu];
        vector<bool>& is_avir_mu = is_avir[mu];
        double***  F_ov_mu    = F_ov[mu];
        double***  F_vv_mu    = F_vv[mu];
        double***  T1_ov_mu   = T1_ov[mu];
        double***  V_oovv_mu  = V_oovv[mu];
        double***  T2_oovv_mu = T2_oovv[mu];

        // Compute P(k/ij) Z_ijk^abc
        Z[mu][ijk.sym]->multiply(T2_ij_a_b.get_block_matrix(ij_abs,mu),
            V_K_bC_e.get_block_matrix(k_abs),0.0,1.0);

        Z[mu][ijk.sym]->multiply(T2_iJ_a_B.get_block_matrix(jk_abs,mu),
            V_k_bC_E.get_block_matrix(i_abs),1.0,-1.0);
        Z[mu][ijk.sym]->multiply(T2_iJ_a_B.get_block_matrix(ik_abs,mu),
            V_k_bC_E.get_block_matrix(j_abs),1.0,1.0);

        Z[mu][ijk.sym]->multiply(V_jk_c_m.get_block_matrix(ij_abs,mu),
            T2_J_aB_i.get_block_matrix(k_abs),1.0,1.0);

        Z[mu][ijk.sym]->multiply(V_jK_c_M.get_block_matrix(jk_abs,mu),
            T2_i_aB_J.get_block_matrix(i_abs),1.0,1.0);
        Z[mu][ijk.sym]->multiply(V_jK_c_M.get_block_matrix(ik_abs,mu),
            T2_i_aB_J.get_block_matrix(j_abs),1.0,-1.0);

        // Compute P(ab) Z_ijk^abc
        W[mu][ijk.sym]->a_b_permutation_1_2(Z[mu][ijk.sym],vvv,v,vv);

        Z[mu][ijk.sym]->multiply(T2_iJ_B_a.get_block_matrix(ik_abs,mu),
            V_k_bc_e.get_block_matrix(j_abs),0.0,1.0);
        Z[mu][ijk.sym]->multiply(T2_iJ_B_a.get_block_matrix(jk_abs,mu),
            V_k_bc_e.get_block_matrix(i_abs),1.0,-1.0);

//        Z[mu][ijk.sym]->multiply(T2_iJ_B_a.get_block_matrix(ik_abs,mu),
//            V_k_bc_e.get_block_matrix(j_abs),1.0,-1.0);
//        Z[mu][ijk.sym]->multiply(T2_iJ_B_a.get_block_matrix(jk_abs,mu),
//            V_k_bc_e.get_block_matrix(i_abs),1.0,1.0);

        W[mu][ijk.sym]->add_c_ab_permutation_1_2(Z[mu][ijk.sym],vvv,v,vv);

        fprintf(outfile,"\n  ijk = %3d  norm = %20.15f",ijk.abs,sqrt(W[mu][ijk.sym]->norm()));


//        // Compute T
//        abc.reset();
//        abc.set_irrep(ijk.sym);
//        while(++abc){
//          int a_sym    = v->get_tuple_irrep(abc.ind_abs[0]);
//          int b_sym    = v->get_tuple_irrep(abc.ind_abs[1]);
//          int c_sym    = v->get_tuple_irrep(abc.ind_abs[2]);
//
//          size_t a_abs = v->get_tuple_abs_index(abc.ind_abs[0]);
//          size_t b_abs = v->get_tuple_abs_index(abc.ind_abs[1]);
//          size_t c_abs = v->get_tuple_abs_index(abc.ind_abs[2]);
//
//          size_t a_rel = v->get_tuple_rel_index(abc.ind_abs[0]);
//          size_t b_rel = v->get_tuple_rel_index(abc.ind_abs[1]);
//          size_t c_rel = v->get_tuple_rel_index(abc.ind_abs[2]);
//
//          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);
//
//          if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_avir_mu[c_abs]){
//            double D_abc = F_vv_mu[a_sym][a_rel][a_rel]
//                         + F_vv_mu[b_sym][b_rel][b_rel]
//                         + F_vv_mu[c_sym][c_rel][c_rel];
//
//            T[mu][ijk.sym]->set(a_sym,a_rel,bc_rel,W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel)/(D_ijk - D_abc));
//          }
//        }

//        // Compute the Energy
//        abc.reset();
//        abc.set_irrep(ijk.sym);
//        while(++abc){
//          int a_sym     = v->get_tuple_irrep(abc.ind_abs[0]);
//
//          size_t a_abs  = v->get_tuple_abs_index(abc.ind_abs[0]);
//          size_t b_abs  = v->get_tuple_abs_index(abc.ind_abs[1]);
//          size_t c_abs  = v->get_tuple_abs_index(abc.ind_abs[2]);
//
//          size_t a_rel  = v->get_tuple_rel_index(abc.ind_abs[0]);
//
//          int bc_sym = vv->get_tuple_irrep(abc.ind_abs[1],abc.ind_abs[2]);
//          size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs[1],abc.ind_abs[2]);
//
//          if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_avir_mu[c_abs]){
//            // E4T
//            E4T[mu]  += W[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) *  T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel) / 36.0;
//
//            // E4ST
//            if((i_sym == a_sym) & (jk_sym == bc_sym)){
//              E4ST[mu] += 0.25 * T1_ov_mu[i_sym][i_rel][a_rel]
//                        * V_oovv_mu[jk_sym][jk_rel][bc_rel]
//                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
//            }
//            // E4DT
//            if((i_sym == a_sym) & (jk_sym == bc_sym)){
//              E4DT[mu] += 0.25 * F_ov_mu[i_sym][i_rel][a_rel]
//                        * T2_oovv_mu[jk_sym][jk_rel][bc_rel]
//                        * T[mu][ijk.sym]->get(a_sym,a_rel,bc_rel);
//            }
//          }
//        }
      }
    }
  }

  double E4T_sum = 0.0;
  double E4ST_sum = 0.0;
  double E4DT_sum = 0.0;
  for(int mu = 0; mu < nurefs; ++mu){
    fprintf(outfile,"\n  E4T(%d)  = %20.15f",mu,E4T[mu]);
    fprintf(outfile,"\n  E4ST(%d) = %20.15f",mu,E4ST[mu]);
    fprintf(outfile,"\n  E4DT(%d) = %20.15f",mu,E4DT[mu]);
    E4T_sum += E4T[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4ST_sum += E4ST[mu] * left_eigenvector[mu] * right_eigenvector[mu];
    E4DT_sum += E4DT[mu] * left_eigenvector[mu] * right_eigenvector[mu];
  }

  // Deallocate Z
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete Z[mu][h];
    }
  }
  release2(Z);

  // Deallocate W
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete W[mu][h];
    }
  }
  release2(W);

  // Deallocate T
  for(int mu = 0; mu < nurefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete T[mu][h];
    }
  }
  release2(T);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



}}  /* End Namespaces */
