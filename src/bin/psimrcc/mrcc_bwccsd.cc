/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "moinfo.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"
#include "utilities.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::update_amps_bwccsd_wrapper()
{
  ptr->update_amps_bwccsd();
}

void CCMRCC::update_amps_bwccsd()
{
  update_t1_amps_bwccsd();
  update_t2_amps_bwccsd();

  zero_internal_delta_amps();

  blas->solve("||Delta_t1||{u}  = t1_delta[o][v]{u} . t1_delta[o][v]{u}");
  blas->solve("||Delta_t1||{u} += t1_delta[O][V]{u} . t1_delta[O][V]{u}");

  blas->solve("||Delta_t2||{u}  = t2_delta[oo][vv]{u} . t2_delta[oo][vv]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[oO][vV]{u} . t2_delta[oO][vV]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[OO][VV]{u} . t2_delta[OO][VV]{u}");

  // Compute the T-AMPS difference
  delta_t1_amps=0.0;
  delta_t2_amps=0.0;
  for(int n=0;n<moinfo->get_ref_size("a");n++){
    delta_t1_amps+=blas->get_scalar("||Delta_t1||",moinfo->get_ref_number("a",n));
    delta_t2_amps+=blas->get_scalar("||Delta_t2||",moinfo->get_ref_number("a",n));
  }
  delta_t1_amps=pow(delta_t1_amps,0.5)/((double)moinfo->get_nrefs());
  delta_t2_amps=pow(delta_t2_amps,0.5)/((double)moinfo->get_nrefs());
}

void CCMRCC::update_t1_amps_bwccsd()
{
  blas->solve("d'1[o][v]{u}  = d1[o][v]{u}");
  blas->solve("d'1[O][V]{u}  = d1[O][V]{u}");

  for(int n=0;n<moinfo->get_ref_size("u");n++){
    double shift = current_energy-Heff[moinfo->get_ref_number("u",n)][moinfo->get_ref_number("u",n)];
    string n_str = to_string(moinfo->get_ref_number("u",n));
    blas->solve("d'1[o][v]{" + n_str + "} += " + to_string(shift));
    blas->solve("d'1[O][V]{" + n_str + "} += " + to_string(shift));
  }

  blas->solve("t1_delta[o][v]{u}  =   t1_eqns[o][v]{u} / d'1[o][v]{u} - t1[o][v]{u}");
  blas->solve("t1_delta[O][V]{u}  =   t1_eqns[O][V]{u} / d'1[O][V]{u} - t1[O][V]{u}");

  blas->solve("t1[o][v]{u} = t1_eqns[o][v]{u} / d'1[o][v]{u}");
  blas->solve("t1[O][V]{u} = t1_eqns[O][V]{u} / d'1[O][V]{u}");

  blas->solve("t1_norm{u}  = t1[o][v]{u} . t1[o][v]{u}");
  blas->solve("t1_norm{u} += t1[O][V]{u} . t1[O][V]{u}");

  zero_t1_internal_amps();

  DEBUGGING(3,
    blas->print("t1[o][v]{u}");
    blas->print("t1[O][V]{u}");
  );
}

void CCMRCC::update_t2_amps_bwccsd()
{
  blas->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
  blas->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
  blas->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");

  for(int n=0;n<moinfo->get_ref_size("u");n++){
    string shift = to_string(current_energy-Heff[moinfo->get_ref_number("u",n)][moinfo->get_ref_number("u",n)]);
    string n_str = to_string(moinfo->get_ref_number("u",n));
    blas->solve("d'2[oo][vv]{" + n_str + "} += " + shift);
    blas->solve("d'2[oO][vV]{" + n_str + "} += " + shift);
    blas->solve("d'2[OO][VV]{" + n_str + "} += " + shift);
  }



  // (a) Compute eq. (20) of J. Chem. Phys. 110, 10275 (1999)
  // Comment : Look at eq. (21) of J. Chem. Phys. 110, 10275 (1999)
  blas->solve("t1_eqns[o][v]{u} += - d1[o][v]{u} * t1[o][v]{u}");
  blas->solve("t1_eqns[O][V]{u} += - d1[O][V]{u} * t1[O][V]{u}");

  // aaaa case
  // (b) Add PijPab (term from a) to the T2 equations 
  blas->solve("t2_eqns[oo][vv]{u} += #1324#   t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2314# - t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #1423# - t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2413#   t1[o][v]{u} X t1_eqns[o][v]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int n=0;n<moinfo->get_ref_size("u");n++){
    int ref = moinfo->get_ref_number("u",n);
    string neg_shift = to_string(-current_energy+Heff[ref][ref]);
    string shift     = to_string(current_energy-Heff[ref][ref]);
    string n_str     = to_string(ref);
    blas->solve("t2_eqns[oo][vv]{" + n_str + "} += #1324# " + neg_shift + "  t1[o][v]{" + n_str + "} X t1[o][v]{" + n_str + "}");
    blas->solve("t2_eqns[oo][vv]{" + n_str + "} += #2314# " + shift  + "  t1[o][v]{" + n_str + "} X t1[o][v]{" + n_str + "}");
  }

  // abab case
  // (b) Add PijPab (term from a) to the T2 equations 
  blas->solve("t2_eqns[oO][vV]{u} += #1324# t1[o][v]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[oO][vV]{u} += #2413# t1[O][V]{u} X t1_eqns[o][v]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int n=0;n<moinfo->get_ref_size("u");n++){
    int ref = moinfo->get_ref_number("u",n);
    string n_str = to_string(ref);
    string neg_shift = to_string(-current_energy+Heff[ref][ref]);
    blas->solve("t2_eqns[oO][vV]{" + n_str + "} += #1324# " + neg_shift + "  t1[o][v]{" + n_str + "} X t1[O][V]{" + n_str + "}");
  }

  // bbbb case
  // (b) Add PijPab (term from a) to the T2 equations 
  blas->solve("t2_eqns[OO][VV]{u} += #1324#   t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #2314# - t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #1423# - t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #2413#   t1[O][V]{u} X t1_eqns[O][V]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int n=0;n<moinfo->get_nunique();n++){
    int ref = moinfo->get_ref_number("u",n);
    string neg_shift = to_string(-current_energy+Heff[ref][ref]);
    string shift     = to_string(current_energy-Heff[ref][ref]);
    string n_str     = to_string(ref);
    blas->solve("t2_eqns[OO][VV]{" + n_str + "} += #1324# " + neg_shift + "  t1[O][V]{" + n_str + "} X t1[O][V]{" + n_str + "}");
    blas->solve("t2_eqns[OO][VV]{" + n_str + "} += #2314# " + shift     + "  t1[O][V]{" + n_str + "} X t1[O][V]{" + n_str + "}");
  }

  blas->solve("t2_delta[oo][vv]{u} = t2_eqns[oo][vv]{u} / d'2[oo][vv]{u} - t2[oo][vv]{u}");
  blas->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u} / d'2[oO][vV]{u} - t2[oO][vV]{u}");
  blas->solve("t2_delta[OO][VV]{u} = t2_eqns[OO][VV]{u} / d'2[OO][VV]{u} - t2[OO][VV]{u}");

  blas->solve("t2[oo][vv]{u} = t2_eqns[oo][vv]{u} / d'2[oo][vv]{u}");
  blas->solve("t2[oO][vV]{u} = t2_eqns[oO][vV]{u} / d'2[oO][vV]{u}");
  blas->solve("t2[OO][VV]{u} = t2_eqns[OO][VV]{u} / d'2[OO][VV]{u}");

  DEBUGGING(3,
    blas->print("t2_eqns[oo][vv]{u}");
    blas->print("t2[oo][vv]{u}");
    blas->print("t2_eqns[oO][vV]{u}");
    blas->print("t2[oO][vV]{u}");
    blas->print("t2_eqns[OO][VV]{u}");
    blas->print("t2[OO][VV]{u}");
  );
}

}} /* End Namespaces */