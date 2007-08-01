/*! \file globals.h
    \ingroup (CCENERGY)
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#ifdef DMALLOC
#include <dmalloc.h>
#endif

/* #define TIME_CCENERGY */

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;

