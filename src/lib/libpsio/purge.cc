/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <map>
#include <sstream>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/workaround.hpp>

using namespace psi;

void PSIO::purge() {
  // for each unit
  for(int u=0; u<psio_unit.size(); ++u) {
    /* Get the number of volumes */
    const int nvol = this->get_numvols(u);

    // get the basename
    char *basename;
    this->get_filename(u, &basename);

    for (int i=0; i<nvol; i++) {
      char* vpath;
      this->get_volpath(u, i, &vpath);

      // get the full filename
      std::ostringstream oss;
      oss << vpath << basename << "." << u;
      const char* fname = oss.str().c_str();

      // unlink
      int errcod = unlink(fname);
      if (errcod != 0) {
        switch (errno) {
          case ENOENT: // no such file? move on
            break;

          default: // else announce to stderr
          {
            std::ostringstream oss;
            oss << "PSIO::purge() -- file=" << fname;
            perror(oss.str().c_str());
          }
        }
      }

      free(vpath);
    }

    free(basename);
  }
}

extern "C" {
  /*!
   ** PSIO_PURGE(): removes all files managed by the default Psio object.
   **
   ** \ingroup PSIO
   */
  int psio_purge() {
    _default_psio_lib_->purge();
    return 1;
  }
}

