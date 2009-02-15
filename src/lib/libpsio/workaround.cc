/*
 * workaround.cc
 *
 *  Created on: Dec 17, 2008
 *      Author: evaleev
 */

#include <libpsio/workaround.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

namespace psi { namespace detail {

  int _OpeN_creatrdwd(const char *path) {
    return ::open(path,O_CREAT|O_RDWR,0644);
  }
  int _OpeN_creatrdwdtrunc(const char *path) {
    return ::open(path,O_CREAT|O_RDWR|O_TRUNC,0644);
  }

}}
