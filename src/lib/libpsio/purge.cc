/*!
 \file
 \ingroup PSIO
 */

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <map>
#include <sstream>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/workaround.hpp>

#define DEBUG 0

using namespace psi;

void PSIO::purge(bool all) {
  // for each unit
  for(int u=0; u<psio_unit.size(); ++u) {
    /* Get the number of volumes */
    const int nvol = this->get_numvols(u);

    // get the basename
    char *basename;
    this->get_filename(u, &basename);

    // purge this file?
    if (all == false && u == PSIF_CHKPT) continue;  // if partial clean requested, do not purge checkpoint file
    bool purge_this = true;
    if (all == false) { // if partial clean, do not purge files that have all volumes residing in current (./) directory
      purge_this = false;
      for (int i=0; i<nvol; i++) {
        char* vpath;
        this->get_volpath(u, i, &vpath);
        if (vpath[0] != '.')
          purge_this = true;
        free(vpath);
      }
    }

    if (purge_this) {
      for (int i = 0; i < nvol; i++) {
        char* vpath;
        this->get_volpath(u, i, &vpath);

        // get the full filename
        char* fname;
        // to avoid bugs in this statement: fname = oss.str().c_str();
        {
          std::ostringstream oss;
          oss << vpath << basename << "." << u;
          std::string fname_str = oss.str();
          const size_t n = fname_str.size();
          fname = new char[n + 1];
          fname_str.copy(fname, n);
          fname[n] = '\0';
        }

        // unlink
        const int errcod = unlink(fname);
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
#if DEBUG
        else {
          std::cerr << "PSIO::purge() -- unlinked " << fname << std::endl;
        }
#endif

        delete[] fname;
        free(vpath);
      }
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

