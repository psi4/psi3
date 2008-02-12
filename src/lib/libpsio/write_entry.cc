/*!
 \file write_entry.cc
 \ingroup PSIO
 */

#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

void PSIO::write_entry(unsigned int unit, char *key, char *buffer, ULI size) {
  psio_address end = PSIO_ZERO;
  write(unit, key, buffer, size, PSIO_ZERO, &end);
}

extern "C" {
  
  /*!
   ** PSIO_WRITE_ENTRY()
   **
   ** \ingroup PSIO
   */
  int psio_write_entry(unsigned int unit, char *key, char *buffer, ULI size) {
    _default_psio_lib_->write_entry(unit, key, buffer, size);
    return 1;
  }
}

