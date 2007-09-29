/*! \file 
    \ingroup (CCHBAR)
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>

namespace psi { namespace cchbar {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

extern "C" {
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
}

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

}} // namespace psi::cchbar
