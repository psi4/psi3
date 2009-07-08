#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "index_iterator.h"

namespace psi{ namespace psimrcc{

using namespace std;

int CCIndexIterator::nirreps_=-1;

CCIndexIterator::CCIndexIterator(string str,int select_irrep)
{
  ccindex_ = blas->get_index(str);
  select = select_irrep;
  if(nirreps_<0) nirreps_ = moinfo->get_nirreps();
  nelements_ = ccindex_->get_nelements();
  element_irrep = ccindex_->get_element_irrep();
  tuples    = ccindex_->get_tuples();
  ind_abs   = new short[nelements_];
  ind_sym   = new int[nelements_];
  reset();
}

CCIndexIterator::CCIndexIterator(CCIndex* index_,int select_irrep)
{
  ccindex_ = index_;
  select = select_irrep;
  if(nirreps_<0) nirreps_ = moinfo->get_nirreps();
  nelements_ = ccindex_->get_nelements();
  element_irrep = ccindex_->get_element_irrep();
  tuples    = ccindex_->get_tuples();
  ind_abs   = new short[nelements_];
  ind_sym   = new int[nelements_];
  reset();
}

CCIndexIterator::~CCIndexIterator()
{
  delete[] ind_abs;
  delete[] ind_sym;
}

void CCIndexIterator::reset()
{
  sym = -1;
  rel = -1;
  abs = -1;
  for(int n = 0; n < nelements_; ++n)
    ind_abs[n] = 0;
  max_index_in_irrep_ = 0;
}

int CCIndexIterator::next_non_empty_irrep(int n)
{
  // special case of one irrep
  if(select != -1){
    if(n < select)
      return(select);
    else
      return(-1);
  }
  for(int h = n + 1; h < nirreps_; ++h)
    if (ccindex_->get_pairpi(h) > 0) return h;
  return -1;
}

bool CCIndexIterator::operator++()
{
  if(rel < max_index_in_irrep_ - 1){
    // Increment index
    ++rel;
    ++abs;
    if(nelements_>0){
      ind_abs[0] = tuples[abs][0];
      ind_sym[0] = element_irrep[0][ind_abs[0]];
    }
    if(nelements_>1){
      ind_abs[1] = tuples[abs][1];
      ind_sym[1] = element_irrep[1][ind_abs[1]];
    }
    if(nelements_>2){
      ind_abs[2] = tuples[abs][2];
      ind_sym[2] = element_irrep[2][ind_abs[2]];
    }
    return(true);
  }else{
    // Check if it reached the last non empty irrep
    sym = next_non_empty_irrep(sym);
    if(sym == -1)
      return(false);
    else{
      max_index_in_irrep_ = ccindex_->get_pairpi(sym);
      rel = 0;
      abs = ccindex_->get_first(sym);
      if(nelements_>0){
        ind_abs[0] = tuples[abs][0];
        ind_sym[0] = element_irrep[0][ind_abs[0]];
      }
      if(nelements_>1){
        ind_abs[1] = tuples[abs][1];
        ind_sym[1] = element_irrep[1][ind_abs[1]];
      }
      if(nelements_>2){
        ind_abs[2] = tuples[abs][2];
        ind_sym[2] = element_irrep[2][ind_abs[2]];
      }
      return(true);
    }
  }
}

bool CCIndexIterator::next_element_in_irrep()
{
  if(rel < max_index_in_irrep_ - 1){
    // Increment index
    ++rel;
    ++abs;
    if(nelements_>0){
      ind_abs[0] = tuples[abs][0];
      ind_sym[0] = element_irrep[0][ind_abs[0]];
    }
    if(nelements_>1){
      ind_abs[1] = tuples[abs][1];
      ind_sym[1] = element_irrep[1][ind_abs[1]];
    }
    if(nelements_>2){
      ind_abs[2] = tuples[abs][2];
      ind_sym[2] = element_irrep[2][ind_abs[2]];
    }
    return(true);
  }else{
    return(false);
  }
}

bool CCIndexIterator::next_irrep()
{
  // Check if it reached the last non empty irrep
  sym = next_non_empty_irrep(sym);
  if(sym == -1)
    return(false);
  else{
    max_index_in_irrep_ = ccindex_->get_pairpi(sym);
    rel = 0;
    abs = ccindex_->get_first(sym);
    if(nelements_>0){
      ind_abs[0] = tuples[abs][0];
      ind_sym[0] = element_irrep[0][ind_abs[0]];
    }
    if(nelements_>1){
      ind_abs[1] = tuples[abs][1];
      ind_sym[1] = element_irrep[1][ind_abs[1]];
    }
    if(nelements_>2){
      ind_abs[2] = tuples[abs][2];
      ind_sym[2] = element_irrep[2][ind_abs[2]];
    }
    return(true);
  }
}

}}
