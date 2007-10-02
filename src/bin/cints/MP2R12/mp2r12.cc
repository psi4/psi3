/*! \file mp2r12.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"moinfo.h"
#include"moinfo_corr.h"
#include"rmp2r12_energy.h"

namespace psi { namespace CINTS {

void mp2r12()
{
  
  init_moinfo();
  init_moinfo_corr();
  switch(UserOptions.reftype) {
  case rhf:
      rmp2r12_energy();
      break;
  default:
      throw std::domain_error("MP2-R12/A energy with specified REFERENCE not implemented");
  }
  cleanup_moinfo_corr();
  cleanup_moinfo();

  return;
}
};};
