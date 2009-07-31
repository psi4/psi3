#include "blas.h"
#include "index_iterator.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

namespace psi{ namespace psimrcc{

double MRCCSD_T::compute_A_ooo_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3)
{
  double value = 0.0;
  int   i_sym  = o->get_tuple_irrep(i_abs);
  int   j_sym  = o->get_tuple_irrep(j_abs);
  int   k_sym  = o->get_tuple_irrep(k_abs);
  int ijk_sym  = i_sym xor j_sym xor k_sym;

  int   x_sym  = v->get_tuple_irrep(x_abs);
  int  jk_sym  = oo->get_tuple_irrep(j_abs,k_abs);

  size_t  x_rel = v->get_tuple_rel_index(x_abs);
  size_t jk_rel = oo->get_tuple_rel_index(j_abs,k_abs);

  if(i_abs == u_abs){
    CCIndexIterator  ef("[vv]",ijk_sym xor x_sym);
    while(++ef){
      int    ef_sym  = vv->get_tuple_irrep(ef.ind_abs[0],ef.ind_abs[1]);
      size_t ef_rel  = vv->get_tuple_rel_index(ef.ind_abs[0],ef.ind_abs[1]);
      if(jk_sym == ef_sym){
        value += 0.25 * T3->get(x_sym,x_rel,ef_rel) * V_oovv[jk_sym][jk_rel][ef_rel];
      }
    }
  }
  return value;
}

double MRCCSD_T::compute_A_ooO_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3)
{
  double value = 0.0;
  int   i_sym  = o->get_tuple_irrep(i_abs);
  int   j_sym  = o->get_tuple_irrep(j_abs);
  int   k_sym  = o->get_tuple_irrep(k_abs);
  int ijk_sym  = i_sym xor j_sym xor k_sym;

  int   x_sym  = v->get_tuple_irrep(x_abs);
  int  jk_sym  = oo->get_tuple_irrep(j_abs,k_abs);

  size_t  x_rel = v->get_tuple_rel_index(x_abs);
  size_t jk_rel = oo->get_tuple_rel_index(j_abs,k_abs);

  if(i_abs == u_abs){
    CCIndexIterator  ef("[vv]",ijk_sym xor x_sym);
    while(++ef){
      int    ef_sym  = vv->get_tuple_irrep(ef.ind_abs[0],ef.ind_abs[1]);
      size_t ef_rel  = vv->get_tuple_rel_index(ef.ind_abs[0],ef.ind_abs[1]);
      if(jk_sym == ef_sym){
        value += T3->get(x_sym,x_rel,ef_rel) * V_oOvV[jk_sym][jk_rel][ef_rel];
      }
    }
  }
  return value;
}

double MRCCSD_T::compute_A_oOO_contribution_to_Heff(int u_abs,int x_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3)
{
  double value = 0.0;
  int   i_sym  = o->get_tuple_irrep(i_abs);
  int   j_sym  = o->get_tuple_irrep(j_abs);
  int   k_sym  = o->get_tuple_irrep(k_abs);
  int ijk_sym  = i_sym xor j_sym xor k_sym;

  int   x_sym  = v->get_tuple_irrep(x_abs);
  int  jk_sym  = oo->get_tuple_irrep(j_abs,k_abs);

  size_t  x_rel = v->get_tuple_rel_index(x_abs);
  size_t jk_rel = oo->get_tuple_rel_index(j_abs,k_abs);

  if(i_abs == u_abs){
    CCIndexIterator  ef("[vv]",ijk_sym xor x_sym);
    while(++ef){
      int    ef_sym  = vv->get_tuple_irrep(ef.ind_abs[0],ef.ind_abs[1]);
      size_t ef_rel  = vv->get_tuple_rel_index(ef.ind_abs[0],ef.ind_abs[1]);
      if(jk_sym == ef_sym){
        value += 0.25 * T3->get(x_sym,x_rel,ef_rel) * V_oovv[jk_sym][jk_rel][ef_rel];
      }
    }
  }
  return value;
}

}} /* End Namespaces */
