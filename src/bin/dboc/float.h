/*! \file float.h
    \ingroup (DBOC)
    \brief Enter brief description of file here 
*/

#include "defines.h"

#ifndef _psi3_dboc_float_h_
#define _psi3_dboc_float_h_

#if LONG_DOUBLE
typedef long double FLOAT;
#else
typedef double FLOAT;
#endif

#if LONG_DOUBLE
#  ifdef AIX
#    define FABS fabsl
#  else
#    define FABS fabs
#  endif
#else
#  define FABS fabs
#endif

#endif
