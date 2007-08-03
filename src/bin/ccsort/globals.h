/*! \file 
    \ingroup (CCSORT)
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;
