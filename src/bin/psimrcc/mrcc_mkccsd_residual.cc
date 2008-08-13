/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include <libmoinfo/libmoinfo.h>
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"
#include <libutil/libutil.h>

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::update_amps_mkccsd_residual_wrapper()
{
  ptr->update_amps_mkccsd_residual();
}

void CCMRCC::update_amps_mkccsd_residual()
{
  zero_internal_amps();
  zero_internal_delta_amps();
  form_similarity_transformed_hamiltonian();
  update_t1_t2_amps_mkccsd_residual();
//   update_t1_amps_mkccsd_residual();
//   update_t2_amps_mkccsd_residual();

  zero_internal_delta_amps();

  blas->solve("||Delta_t1||{u}  = t1_delta[o][v]{u} . t1_delta[o][v]{u}");
  blas->solve("||Delta_t1||{u} += t1_delta[O][V]{u} . t1_delta[O][V]{u}");

  blas->solve("||Delta_t2||{u}  = t2_delta[oo][vv]{u} . t2_delta[oo][vv]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[oO][vV]{u} . t2_delta[oO][vV]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[OO][VV]{u} . t2_delta[OO][VV]{u}");

  // Compute the T-AMPS difference
  delta_t1_amps=0.0;
  delta_t2_amps=0.0;
  for(int n=0;n<moinfo->get_nunique();n++){
    int m = moinfo->get_ref_number("u",n);
    delta_t1_amps+=blas->get_scalar("||Delta_t1||",m);
    delta_t2_amps+=blas->get_scalar("||Delta_t2||",m);
  }
  delta_t1_amps=pow(delta_t1_amps,0.5)/((double)moinfo->get_nunique());
  delta_t2_amps=pow(delta_t2_amps,0.5)/((double)moinfo->get_nunique());
}

void CCMRCC::form_similarity_transformed_hamiltonian()
{
  blas->solve("t1_eqns[o][v]{u}   += - d1[o][v]{u} * t1[o][v]{u}");
  blas->solve("t1_eqns[O][V]{u}   += - d1[O][V]{u} * t1[O][V]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += - d2[oo][vv]{u} * t2[oo][vv]{u}");
  blas->solve("t2_eqns[oO][vV]{u} += - d2[oO][vV]{u} * t2[oO][vV]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += - d2[OO][VV]{u} * t2[OO][VV]{u}");
  zero_internal_delta_amps();
//   for(int i=0;i<moinfo->get_nunique();i++){
//     int unique_i = moinfo->get_ref_number("u",i);
//     string i_str = to_string(unique_i);
//     string factor = to_string(fabs(eigenvector[unique_i]));
//     blas->solve("Mk1[o][v]{" + i_str + "} = " + factor + " t1_eqns[o][v]{" + i_str + "}");
//     blas->solve("Mk1[O][V]{" + i_str + "} = " + factor + " t1_eqns[O][V]{" + i_str + "}");
//     blas->solve("Mk2[oo][vv]{" + i_str + "} = " + factor + " t2_eqns[oo][vv]{" + i_str + "}");
//     blas->solve("Mk2[oO][vV]{" + i_str + "} = " + factor + " t2_eqns[oO][vV]{" + i_str + "}");
//     blas->solve("Mk2[OO][VV]{" + i_str + "} = " + factor + " t2_eqns[OO][VV]{" + i_str + "}");
//   }
//   blas->solve("t1_eqns[o][v]{u} = Mk1[o][v]{u}");
//   blas->solve("t1_eqns[O][V]{u} = Mk1[O][V]{u}");
//   blas->solve("t2_eqns[oo][vv]{u} = Mk2[oo][vv]{u}");
//   blas->solve("t2_eqns[oO][vV]{u} = Mk2[oO][vV]{u}");
//   blas->solve("t2_eqns[OO][VV]{u} = Mk2[OO][VV]{u}");
}

void CCMRCC::update_t1_t2_amps_mkccsd_residual()
{
  blas->solve("d'1[o][v]{u}  = d1[o][v]{u}");
  blas->solve("d'1[O][V]{u}  = d1[O][V]{u}");

  blas->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
  blas->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
  blas->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");

  for(int n=0;n<moinfo->get_nunique();n++){
    int m = moinfo->get_ref_number("u",n);
    string shift = to_string(current_energy-Heff[m][m]);
    blas->solve("d'1[o][v]{" + to_string(m) + "} += " + shift);
    blas->solve("d'1[O][V]{" + to_string(m) + "} += " + shift);
    blas->solve("d'2[oo][vv]{" + to_string(m) + "} += " + shift);
    blas->solve("d'2[oO][vV]{" + to_string(m) + "} += " + shift);
    blas->solve("d'2[OO][VV]{" + to_string(m) + "} += " + shift);
  }

  for(int i=0;i<moinfo->get_nunique();i++){
    int unique_i = moinfo->get_ref_number("u",i);
    string i_str = to_string(unique_i);
    // Form the coupling terms
    for(int j=0;j<moinfo->get_nrefs();j++){
      int unique_j = moinfo->get_ref_number("a",j);
      string j_str = to_string(unique_j);
//       if(fabs(eigenvector[unique_i])>0.0){
        blas->set_scalar("factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_i]);
        blas->set_scalar("neg_factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_i]);
//       }else{
//         CCMatrix::set_scalar("factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]);
//         CCMatrix::set_scalar("neg_factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]);
//       }
      if(unique_i!=j){
        if(j==unique_j){
          blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
          blas->solve("t1_eqns[o][v]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[o][v]{" + i_str + "}");

          blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
          blas->solve("t1_eqns[O][V]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[O][V]{" + i_str + "}");
        }else{
          blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
          blas->solve("t1_eqns[o][v]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[o][v]{" + i_str + "}");

          blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
          blas->solve("t1_eqns[O][V]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[O][V]{" + i_str + "}");
        }
      }
    }

    // Update t1 for reference i
    blas->solve("t1_delta[o][v]{" + i_str + "}  =   t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "} - t1[o][v]{" + i_str + "}");
    blas->solve("t1_delta[O][V]{" + i_str + "}  =   t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "} - t1[O][V]{" + i_str + "}");

    blas->solve("t1[o][v]{" + i_str + "} += t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "}");
    blas->solve("t1[O][V]{" + i_str + "} += t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "}");
    zero_internal_amps();

    // Add the contribution from the other references
    for(int j=0;j<moinfo->get_nrefs();j++){
      int unique_j = moinfo->get_ref_number("a",j);
      string j_str = to_string(unique_j);
      blas->set_scalar("factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_i]);
      //       if(fabs(eigenvector[unique_i])>0.0)
//         CCMatrix::set_scalar("factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]);
//       else
//         CCMatrix::set_scalar("factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]);
      if(unique_i!=j){
        if(j==unique_j){
          // aaaa case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");

          // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");

          // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");

          // P(ij)t_i^a(mu)t_j^b(mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += - t2[oo][vv]{" + i_str + "}");
  
          blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");
  
          // abab case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "}  = t2[oO][vV]{" + j_str + "}");
  
          // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[O][V]{" + j_str + "}");
  
          // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + j_str + "} X t1[O][V]{" + i_str + "}");
  
          // P(iJ)t_i^a(mu)t_J^B(mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += - t2[oO][vV]{" + i_str + "}");

          blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");

          // bbbb case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");

          // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");

          // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");

          // P(ij)t_i^a(mu)t_j^b(mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += - t2[OO][VV]{" + i_str + "}");

          blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
        }else{
          // aaaa case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");

          // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");

          // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");

          // P(ij)t_i^a(mu)t_j^b(mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
          blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[oo][vv]{" + i_str + "} += - t2[oo][vv]{" + i_str + "}");

          blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");
  
          // abab case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "}  = #2143# t2[oO][vV]{" + j_str + "}");
  
          // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[o][v]{" + j_str + "}");
  
          // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[O][V]{" + j_str + "} X t1[O][V]{" + i_str + "}");

          // P(iJ)t_i^a(mu)t_J^B(mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[oO][vV]{" + i_str + "} += - t2[oO][vV]{" + i_str + "}");

          blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");

          // bbbb case
          // + t_ij^ab(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");

          // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");

          // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");

          // P(ij)t_i^a(mu)t_j^b(mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
          blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");

          // - t_ij^ab(mu)
          blas->solve("Mk2[OO][VV]{" + i_str + "} += - t2[OO][VV]{" + i_str + "}");

          blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
        }

      }
    }
    blas->solve("t2_delta[oo][vv]{" + i_str + "} = t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "} - t2[oo][vv]{" + i_str + "}");
    blas->solve("t2_delta[oO][vV]{" + i_str + "} = t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "} - t2[oO][vV]{" + i_str + "}");
    blas->solve("t2_delta[OO][VV]{" + i_str + "} = t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "} - t2[OO][VV]{" + i_str + "}");

    blas->solve("t2[oo][vv]{" + i_str + "} += t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "}");
    blas->solve("t2[oO][vV]{" + i_str + "} += t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "}");
    blas->solve("t2[OO][VV]{" + i_str + "} += t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "}");
    zero_internal_amps();
  }



  blas->solve("t1_norm{u}  = t1[o][v]{u} . t1[o][v]{u}");
  blas->solve("t1_norm{u} += t1[O][V]{u} . t1[O][V]{u}");



  DEBUGGING(3,
    blas->print("t2_eqns[oo][vv]{u}");
    blas->print("t2[oo][vv]{u}");
    blas->print("t2_eqns[oO][vV]{u}");
    blas->print("t2[oO][vV]{u}");
    blas->print("t2_eqns[OO][VV]{u}");
    blas->print("t2[OO][VV]{u}");
  )
}

// void CCMRCC::update_t1_amps_mkccsd_residual()
// {
//   blas->solve("d'1[o][v]{u}  = d1[o][v]{u}");
//   blas->solve("d'1[O][V]{u}  = d1[O][V]{u}");
// 
//   for(int n=0;n<moinfo->get_nunique();n++){
//     int m = moinfo->get_ref_number("u",n);
//     string shift = to_string(moinfo->get_denominator_shift()+current_energy-Heff[m][m]);
//     blas->solve("d'1[o][v]{" + to_string(m) + "} += " + shift);
//     blas->solve("d'1[O][V]{" + to_string(m) + "} += " + shift);
//   }
// 
//   for(int i=0;i<moinfo->get_nunique();i++){
//     int unique_i = moinfo->get_ref_number("u",i);
//     string i_str = to_string(unique_i);
//     // Form the coupling terms
//     for(int j=0;j<moinfo->get_nrefs();j++){
//       int unique_j = moinfo->get_ref_number("a",j);
//       string j_str = to_string(unique_j);
// //       if(eigenvector[unique_i]>0.0){
//         CCMatrix::set_scalar("factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_i]);
//         CCMatrix::set_scalar("neg_factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_i]);
// //       }
// //       else{
// //         CCMatrix::set_scalar("factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]);
// //         CCMatrix::set_scalar("neg_factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]);
// //       }
//       if(unique_i!=j){
//         if(j==unique_j){
//           blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
//           blas->solve("t1_eqns[o][v]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[o][v]{" + i_str + "}");
// 
//           blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
//           blas->solve("t1_eqns[O][V]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[O][V]{" + i_str + "}");
//         }else{
//           blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
//           blas->solve("t1_eqns[o][v]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[o][v]{" + i_str + "}");
// 
//           blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
//           blas->solve("t1_eqns[O][V]{" + i_str + "} += neg_factor_mk{" + j_str + "} t1[O][V]{" + i_str + "}");
//         }
//       }
//     }
//     // Update t1 for reference i
//     blas->solve("t1_delta[o][v]{" + i_str + "}  =   t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "}");
//     blas->solve("t1_delta[O][V]{" + i_str + "}  =   t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "}");
// 
//     blas->solve("t1[o][v]{" + i_str + "} += t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "}");
//     blas->solve("t1[O][V]{" + i_str + "} += t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "}");
//     zero_internal_amps();
//   }
// }
// 
// void CCMRCC::update_t2_amps_mkccsd_residual()
// {
//   blas->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
//   blas->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
//   blas->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");
// 
//   for(int n=0;n<moinfo->get_nunique();n++){
//     int m = moinfo->get_ref_number("u",n);
//     string shift = to_string(moinfo->get_denominator_shift()+current_energy-Heff[m][m]);
//     blas->solve("d'2[oo][vv]{" + to_string(m) + "} += " + shift);
//     blas->solve("d'2[oO][vV]{" + to_string(m) + "} += " + shift);
//     blas->solve("d'2[OO][VV]{" + to_string(m) + "} += " + shift);
//   }
// 
//   for(int i=0;i<moinfo->get_nunique();i++){
//     int unique_i = moinfo->get_ref_number("u",i);
//     string i_str = to_string(unique_i);
//     for(int j=0;j<moinfo->get_nrefs();j++){
//       int unique_j = moinfo->get_ref_number("a",j);
//       string j_str = to_string(unique_j);
// //       if(eigenvector[unique_i]>0.0)
//         CCMatrix::set_scalar("factor_mk",unique_j,Heff[unique_i][j]*eigenvector[j]/eigenvector[unique_j]);
// //       else
// //         CCMatrix::set_scalar("factor_mk",unique_j,-Heff[unique_i][j]*eigenvector[j]);
//       if(unique_i!=j){
//         if(j==unique_j){
//           // aaaa case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
// 
//           // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
// 
//           // P(ij)t_i^a(mu)t_j^b(mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
// 
//           // - t_ij^ab(mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += - t2[oo][vv]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");
// 
//           // abab case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "}  = t2[oO][vV]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[O][V]{" + j_str + "}");
// 
//           // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + j_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // P(iJ)t_i^a(mu)t_J^B(mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // - t_ij^ab(mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += - t2[oO][vV]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");
// 
//           // bbbb case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
// 
//           // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
// 
//           // P(ij)t_i^a(mu)t_j^b(mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // + t_ij^ab(mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += - t2[OO][VV]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
//         }else{
//           // aaaa case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
// 
//           // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
// 
//           // P(ij)t_i^a(mu)t_j^b(mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
// 
//           // - t_ij^ab(mu)
//           blas->solve("Mk2[oo][vv]{" + i_str + "} += - t2[oo][vv]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");
// 
//           // abab case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "}  = #2143# t2[oO][vV]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[o][v]{" + j_str + "}");
// 
//           // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[O][V]{" + j_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // P(iJ)t_i^a(mu)t_J^B(mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // + t_ij^ab(mu)
//           blas->solve("Mk2[oO][vV]{" + i_str + "} += - t2[oO][vV]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");
// 
//           // bbbb case
//           // + t_ij^ab(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");
// 
//           // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
// 
//           // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
// 
//           // P(ij)t_i^a(mu)t_j^b(mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
// 
//           // + t_ij^ab(mu)
//           blas->solve("Mk2[OO][VV]{" + i_str + "} += - t2[OO][VV]{" + i_str + "}");
// 
//           blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
//         }
// 
//       }
//     }
//     blas->solve("t2_delta[oo][vv]{" + i_str + "} = t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "}");
//     blas->solve("t2_delta[oO][vV]{" + i_str + "} = t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "}");
//     blas->solve("t2_delta[OO][VV]{" + i_str + "} = t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "}");
// 
//     blas->solve("t2[oo][vv]{" + i_str + "} += t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "}");
//     blas->solve("t2[oO][vV]{" + i_str + "} += t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "}");
//     blas->solve("t2[OO][VV]{" + i_str + "} += t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "}");
//     zero_internal_amps();
//   }
// 
//   if(moinfo->get_debug()>3){
//     blas->print("t2_eqns[oo][vv]{u}");
//     blas->print("t2[oo][vv]{u}");
//     blas->print("t2_eqns[oO][vV]{u}");
//     blas->print("t2[oO][vV]{u}");
//     blas->print("t2_eqns[OO][VV]{u}");
//     blas->print("t2[OO][VV]{u}");
//   }
// }


}} /* End Namespaces */
