#include <stdio.h>
#include <ccfiles.h>
#include <dpd.h>
#include "MOInfo.h"

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN struct MOInfo moinfo;

EXTERN int *triotot, ***triorb;

