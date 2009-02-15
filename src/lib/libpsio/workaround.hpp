/*
 * workaround.hpp
 *
 *  Created on: Dec 17, 2008
 *      Author: evaleev
 */

#ifndef _psi_src_lib_libpsio_workaround_hpp_
#define _psi_src_lib_libpsio_workaround_hpp_

// wrappers around open() designed to avoid including broken fcntl.h in the file
// containing PSIO::open() (reported to include macro redefining open to open64)
namespace psi { namespace detail {

  int _OpeN_creatrdwd(const char *path);
  int _OpeN_creatrdwdtrunc(const char *path);

}} // end of namespace psi::detail

#endif /* header guard */
