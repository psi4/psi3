/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
extern "C" {
#include <libpsio/psio.h>
}
#include <libpsio/psio.hpp>

using namespace psi;

/* Definition of global data */
PSIO* psi::_default_psio_lib_ = 0;
int PSIO::_error_exit_code_ = 1;
extern "C" {
  psio_address PSIO_ZERO = { 0, 0 };
}

PSIO::PSIO() : psio_unit(PSIO_MAXUNIT) {
  state_ = 1;
  
  for (int i=0; i < psio_unit.size(); i++) {
#ifdef PSIO_STATS
    psio_unit[i].readlen = psio_unit[i].writlen = 0;
#endif      
    psio_unit[i].numvols = 0;
    for (int j=0; j < PSIO_MAXVOL; j++) {
      psio_unit[i].vol[j].path = NULL;
      psio_unit[i].vol[j].stream = -1;
    }
    psio_unit[i].toclen = 0;
    psio_unit[i].toc = NULL;
  }
}

extern "C" {
  
  /*!
   ** PSIO_INIT(): Allocates global memory needed by the I/O routines.
   **
   ** No arguments.
   **
   ** \ingroup PSIO
   */

  int psio_init(void) {
    if (!_default_psio_lib_) {
      _default_psio_lib_ = new PSIO;
    }
    
    return 1;
  }
  
  /*!
   ** PSIO_STATE(): Returns state of the library (1=initialized, 0=noninitialized).
   **
   ** No arguments.
   **
   ** \ingroup PSIO
   */

  int psio_state() {
    return _default_psio_lib_->state();
  }

}
