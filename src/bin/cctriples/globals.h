/*! \file 
    \ingroup (CCTRIPLES)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

EXTERN int *triotot, ***triorb;

